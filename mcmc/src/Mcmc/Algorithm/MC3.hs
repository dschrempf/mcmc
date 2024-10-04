{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Algorithm.MC3
-- Description :  Metropolis-coupled Markov chain Monte Carlo algorithm
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 23 15:20:33 2020.
--
-- The Metropolis-coupled Markov chain Monte Carlo ('MC3') algorithm.
--
-- Also known as parallel tempering.
--
-- Like any other parallel MCMC algorithm, the 'MC3' algorithm is essentially an
-- 'Mcmc.Algorithm.MHG.MHG' algorithm on the product space of all parallel
-- chains.
--
-- For example, see
--
-- - Geyer, C. J., Markov chain monte carlo maximum likelihood, Computing
--   Science and Statistics, Proceedings of the 23rd Symposium on the Interface,
--   (1991).
--
-- - Altekar, G., Dwarkadas, S., Huelsenbeck, J. P., & Ronquist, F., Parallel
--   metropolis coupled markov chain monte carlo for bayesian phylogenetic
--   inference, Bioinformatics, 20(3), 407–415 (2004).
module Mcmc.Algorithm.MC3
  ( -- * Definitions
    NChains (..),
    SwapPeriod (..),
    NSwaps (..),
    MC3Settings (..),
    MHGChains,
    ReciprocalTemperatures,

    -- * Metropolis-coupled Markov chain Monte Carlo algorithm
    MC3 (..),
    mc3,
    mc3Save,
    mc3Load,
  )
where

import Codec.Compression.GZip
import Control.Concurrent.Async hiding (link)
import Control.Monad
import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List
import qualified Data.Map.Strict as M
import Data.Time
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Word
import Mcmc.Acceptance
import Mcmc.Algorithm
import Mcmc.Algorithm.MHG
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Save
import Mcmc.Chain.Trace
import Mcmc.Cycle
import Mcmc.Internal.Random
import Mcmc.Internal.Shuffle
import Mcmc.Likelihood
import Mcmc.Monitor
import Mcmc.Posterior
import Mcmc.Prior
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log hiding (sum)
import System.Random.Stateful
import Text.Printf

-- | Total number of parallel chains.
--
-- Must be two or larger.
newtype NChains = NChains {fromNChains :: Int}
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''NChains)

-- | The period of proposing state swaps between chains.
--
-- Must be one or larger.
newtype SwapPeriod = SwapPeriod {fromSwapPeriod :: Int}
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SwapPeriod)

-- | The number of proposed swaps at each swapping event.
--
-- Must be in @[1, NChains - 1]@.
newtype NSwaps = NSwaps {fromNSwaps :: Int}
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''NSwaps)

-- | MC3 settings.
data MC3Settings = MC3Settings
  { -- | The number of chains has to be larger equal two.
    mc3NChains :: NChains,
    mc3SwapPeriod :: SwapPeriod,
    mc3NSwaps :: NSwaps
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''MC3Settings)

-- | Vector of MHG chains.
type MHGChains a = V.Vector (MHG a)

-- | Vector of reciprocal temperatures.
type ReciprocalTemperatures = U.Vector Double

data SavedMC3 a = SavedMC3
  { savedMC3Settings :: MC3Settings,
    savedMC3Chains :: V.Vector (SavedChain a),
    savedMC3ReciprocalTemperatures :: ReciprocalTemperatures,
    savedMC3Iteration :: Int,
    savedMC3SwapAcceptance :: Acceptances Int,
    savedMC3Generator :: (Word64, Word64)
  }
  deriving (Eq, Show)

$(deriveJSON defaultOptions ''SavedMC3)

toSavedMC3 ::
  MC3 a ->
  IO (SavedMC3 a)
toSavedMC3 (MC3 s mhgs bs i ac g) = do
  scs <- V.mapM (toSavedChain . fromMHG) mhgs
  g' <- saveGen g
  return $ SavedMC3 s scs bs i ac g'

fromSavedMC3 ::
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  SavedMC3 a ->
  IO (MC3 a)
fromSavedMC3 pr lh cc mn (SavedMC3 s scs bs i ac g') = do
  mhgs <-
    V.fromList
      <$> sequence
        [ MHG <$> fromSavedChain pf lf cc mn sc
          | (sc, pf, lf) <- zip3 (V.toList scs) prs lhs
        ]
  g <- loadGen g'
  return $ MC3 s mhgs bs i ac g
  where
    prs = map (heatFunction pr) $ U.toList bs
    lhs = map (heatFunction lh) $ U.toList bs

-- | The MC3 algorithm.
data MC3 a = MC3
  { mc3Settings :: MC3Settings,
    -- | The first chain is the cold chain with temperature 1.0.
    mc3MHGChains :: MHGChains a,
    -- | Vector of reciprocal temperatures.
    mc3ReciprocalTemperatures :: ReciprocalTemperatures,
    -- | Current iteration.
    mc3Iteration :: Int,
    -- | Number of accepted and rejected swaps.
    mc3SwapAcceptances :: Acceptances Int,
    mc3Generator :: IOGenM StdGen
  }

instance (ToJSON a) => Algorithm (MC3 a) where
  aName = const "Metropolis-coupled Markov chain Monte Carlo (MC3)"
  aIteration = mc3Iteration
  aIsInvalidState = mc3IsInvalidState
  aIterate = mc3Iterate
  aAutoTune = mc3AutoTune
  aResetAcceptance = mc3ResetAcceptance
  aCleanAfterBurnIn = mc3CleanAfterBurnIn
  aSummarizeCycle = mc3SummarizeCycle
  aOpenMonitors = mc3OpenMonitors
  aExecuteMonitors = mc3ExecuteMonitors
  aStdMonitorHeader = mc3StdMonitorHeader
  aCloseMonitors = mc3CloseMonitors
  aSave = mc3Save

heatFunction ::
  -- Cold Function.
  (a -> Log Double) ->
  -- Reciprocal temperature.
  Double ->
  -- The heated prior or likelihood function
  (a -> Log Double)
heatFunction f b
  | b <= 0 = error "heatFunction: Reciprocal temperature is zero or negative."
  | b == 1.0 = f
  | otherwise = (** b') . f
  where
    b' = Exp $ log b

--  The prior and likelihood values of the current link are updated.
--
-- NOTE: The trace is not changed! In particular, the prior and likelihood
-- values are not updated for any link of the trace, and no new link is added to
-- the trace.
setReciprocalTemperature ::
  -- Cold prior function.
  PriorFunction a ->
  -- Cold likelihood function.
  LikelihoodFunction a ->
  -- New reciprocal temperature.
  Double ->
  MHG a ->
  MHG a
setReciprocalTemperature coldPrf coldLhf b a =
  MHG $
    c
      { priorFunction = prf',
        likelihoodFunction = lhf',
        link = Link x (prf' x) (lhf' x)
      }
  where
    c = fromMHG a
    -- We need twice the amount of computations compared to taking the power
    -- after calculating the posterior (pr x * lh x) ** b'. However, I don't
    -- think this is a serious problem.
    --
    -- To minimize computations, it is key to avoid modification of the
    -- reciprocal temperature for the cold chain.
    prf' = heatFunction coldPrf b
    lhf' = heatFunction coldLhf b
    x = state $ link c

initMHG ::
  -- Cold prior function.
  PriorFunction a ->
  -- Cold likelihood function.
  LikelihoodFunction a ->
  -- Index of MHG chain.
  Int ->
  -- Reciprocal temperature.
  Double ->
  MHG a ->
  IO (MHG a)
initMHG prf lhf i beta a
  | i < 0 = error "initMHG: Chain index negative."
  -- Do nothing for the cold chain.
  | i == 0 = return $ MHG c
  | otherwise = do
      -- We have to push the current link in the trace, since it is not set by
      -- 'setReciprocalTemperature'. The other links in the trace are still
      -- pointing to the link of the cold chain, but this has no effect.
      t' <- pushT l t
      return $ MHG $ c {trace = t'}
  where
    a' = setReciprocalTemperature prf lhf beta a
    c = fromMHG a'
    l = link c
    t = trace c

-- | Initialize an MC3 algorithm with a given number of chains.
--
-- Call 'error' if:
--
-- - The number of chains is one or lower.
--
-- - The swap period is zero or negative.
mc3 ::
  MC3Settings ->
  Settings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  InitialState a ->
  StdGen ->
  IO (MC3 a)
mc3 sMc3 s pr lh cc mn i0 g
  | n < 2 = error "mc3: The number of chains must be two or larger."
  | sp < 1 = error "mc3: The swap period must be strictly positive."
  | sn < 1 || sn > n - 1 = error "mc3: The number of swaps must be in [1, NChains - 1]."
  | otherwise = do
      -- Split random number generator.
      let gs = take (n + 1) $ unfoldr (Just . split) g
      -- Prepare MHG chains.
      cs <- V.mapM (mhg s pr lh cc mn i0) (V.fromList $ tail gs)
      hcs <- V.izipWithM (initMHG pr lh) (V.convert bs) cs
      -- Do not reuse the initial generator.
      gm <- newIOGenM $ head gs
      return $ MC3 sMc3 hcs bs 0 (emptyA [0 .. n - 2]) gm
  where
    n = fromNChains $ mc3NChains sMc3
    sp = fromSwapPeriod $ mc3SwapPeriod sMc3
    sn = fromNSwaps $ mc3NSwaps sMc3
    -- NOTE: The initial choice of reciprocal temperatures is based on a few
    -- tests but otherwise pretty arbitrary.
    --
    -- NOTE: Have to 'take n' elements, because vectors are not as lazy as
    -- lists.
    bs = U.fromList $ take n $ iterate (* 0.97) 1.0

mc3Fn :: AnalysisName -> FilePath
mc3Fn (AnalysisName nm) = nm ++ ".mcmc.mc3"

-- | Save an MC3 algorithm.
mc3Save ::
  (ToJSON a) =>
  AnalysisName ->
  MC3 a ->
  IO ()
mc3Save nm a = do
  savedMC3 <- toSavedMC3 a
  BL.writeFile (mc3Fn nm) $ compress $ encode savedMC3

-- | Load an MC3 algorithm.
--
-- Also create a backup of the save.
--
-- See 'Mcmc.Mcmc.mcmcContinue'.
mc3Load ::
  (FromJSON a) =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  AnalysisName ->
  IO (MC3 a)
mc3Load pr lh cc mn nm = do
  savedMC3 <- eitherDecode . decompress <$> BL.readFile fn
  either error (fromSavedMC3 pr lh cc mn) savedMC3
  where
    -- fnBak = mc3Fn $ AnalysisName $ (fromAnalysisName nm ++ ".bak")
    fn = mc3Fn nm

-- I call the chains left and right, because it is easy to think about them as
-- being left and right. Of course, the left chain may also have a larger index
-- than the right chain.
swapWith ::
  -- Index i>=0 of left chain.
  Int ->
  -- Index j>=0, j/=i of right chain.
  Int ->
  MHGChains a ->
  (MHGChains a, Posterior)
swapWith i j xs
  | i < 0 = error "swapWith: Left index is negative."
  | j < 0 = error "swapWith: Right index is negative."
  | i == j = error "swapWith: Indices are equal."
  | otherwise = (xs', q)
  where
    -- Gather information from current chains.
    cl = fromMHG $ xs V.! i
    cr = fromMHG $ xs V.! j
    ll = link cl
    lr = link cr
    prl = prior ll
    prr = prior lr
    lhl = likelihood ll
    lhr = likelihood lr
    -- Swap the states.
    xl' = state lr
    xr' = state ll
    -- Compute new priors and likelihoods.
    prl' = priorFunction cl xl'
    prr' = priorFunction cr xr'
    lhl' = likelihoodFunction cl xl'
    lhr' = likelihoodFunction cr xr'
    -- Set the new links and the proposed state.
    ll' = Link xl' prl' lhl'
    lr' = Link xr' prr' lhr'
    cl' = cl {link = ll'}
    cr' = cr {link = lr'}
    xs' = xs V.// [(i, MHG cl'), (j, MHG cr')]
    -- Compute the Metropolis ratio.
    nominator = prl' * prr' * lhl' * lhr'
    denominator = prl * prr * lhl * lhr
    q = nominator / denominator

mc3ProposeSwap ::
  MC3 a ->
  -- Index of left chain.
  Int ->
  IO (MC3 a)
mc3ProposeSwap a i = do
  let cs = mc3MHGChains a
  -- -- Debug.
  -- prL = prior $ link $ fromMHG $ cs V.! i
  -- prR = prior $ link $ fromMHG $ cs V.! (i+1)
  -- lhL = likelihood $ link $ fromMHG $ cs V.! i
  -- lhR = likelihood $ link $ fromMHG $ cs V.! (i+1)
  -- 1. Sample new state and get the Metropolis ratio.
  let (!y, !r) = swapWith i (i + 1) cs
  -- 2. Accept or reject.
  accept <- mhgAccept r g
  if accept
    then do
      -- -- Debug.
      -- traceIO $ "Swap accepted: " <> show i <> " <-> " <> show (i+1)
      -- let prL' = prior $ link $ fromMHG $ y V.! i
      --     prR' = prior $ link $ fromMHG $ y V.! (i+1)
      --     lhL' = likelihood $ link $ fromMHG $ y V.! i
      --     lhR' = likelihood $ link $ fromMHG $ y V.! (i+1)
      -- traceIO $ "Log priors (left, right, before swap): " <> show (ln prL) <> " " <> show (ln prR)
      -- traceIO $ "Log priors (left, right, after swap): " <> show (ln prL') <> " " <> show (ln prR')
      -- traceIO $ "Log likelihoods (left, right, before swap): " <> show (ln lhL) <> " " <> show (ln lhR)
      -- traceIO $ "Log likelihood (left, right, after swap): " <> show (ln lhL') <> " " <> show (ln lhR')
      let !ac' = pushAccept Nothing i (mc3SwapAcceptances a)
      return $ a {mc3MHGChains = y, mc3SwapAcceptances = ac'}
    else do
      let !ac' = pushReject Nothing i (mc3SwapAcceptances a)
      return $ a {mc3SwapAcceptances = ac'}
  where
    g = mc3Generator a

mc3IsInvalidState :: (ToJSON a) => MC3 a -> Bool
mc3IsInvalidState a = V.any aIsInvalidState mhgs
  where
    mhgs = mc3MHGChains a

-- NOTE: 'mc3Iterate' is actually not parallel, but concurrent because of the IO
-- constraint of the mutable trace.
mc3Iterate ::
  (ToJSON a) =>
  IterationMode ->
  ParallelizationMode ->
  MC3 a ->
  IO (MC3 a)
mc3Iterate m pm a = do
  -- 1. Maybe propose swaps.
  --
  -- NOTE: Swaps have to be proposed first, because the traces are automatically
  -- updated at step 2.
  let s = mc3Settings a
  a' <-
    if mc3Iteration a `mod` fromSwapPeriod (mc3SwapPeriod s) == 0
      then do
        let n = V.length $ mc3MHGChains a
            is = [0 .. n - 2]
            ns = fromNSwaps $ mc3NSwaps s
        is' <- shuffle is $ mc3Generator a
        foldM mc3ProposeSwap a (take ns is')
      else return a
  -- 2. Iterate all chains and increment iteration.
  mhgs <- case pm of
    Sequential -> V.mapM (aIterate m pm) (mc3MHGChains a')
    Parallel ->
      -- Go via a list, and use 'forkIO' ("Control.Concurrent.Async").
      V.fromList <$> mapConcurrently (aIterate m pm) (V.toList $ mc3MHGChains a')
  let i = mc3Iteration a'
  return $ a' {mc3MHGChains = mhgs, mc3Iteration = succ i}

tuneBeta ::
  -- The old reciprocal temperatures are needed to retrieve the old ratios.
  ReciprocalTemperatures ->
  -- Index i of left chain. Change the reciprocal temperature of chain (i+1).
  Int ->
  -- Exponent xi of the reciprocal temperature ratio.
  Double ->
  -- The new reciprocal temperatures are updated incrementally using the
  -- reciprocal temperature ratios during the fold (see 'mc3AutoTune' below).
  ReciprocalTemperatures ->
  ReciprocalTemperatures
tuneBeta bsOld i xi bsNew = bsNew U.// [(j, brNew)]
  where
    j = i + 1
    blOld = bsOld U.! i
    brOld = bsOld U.! j
    blNew = bsNew U.! i
    -- The new ratio is in (0,1).
    rNew = (brOld / blOld) ** xi
    brNew = blNew * rNew

mc3AutoTune :: (ToJSON a) => TuningType -> Int -> MC3 a -> IO (MC3 a)
mc3AutoTune b l a = do
  -- 1. Auto tune all chains.
  mhgs' <- V.mapM (aAutoTune b l) $ mc3MHGChains a
  -- 2. Auto tune temperatures.
  let optimalRate = getOptimalRate PDimensionUnknown
      mCurrentRates = acceptanceRates $ mc3SwapAcceptances a
      -- We assume that the acceptance rate of state swaps between two chains is
      -- roughly proportional to the ratio of the temperatures of the chains.
      -- Hence, we focus on temperature ratios, actually reciprocal temperature
      -- ratios, which is the same. Also, by working with ratios in (0,1) of
      -- neighboring chains, we ensure the monotonicity of the reciprocal
      -- temperatures.
      --
      -- The factor (1/2) was determined by a few tests and is otherwise
      -- absolutely arbitrary.
      xi i = case mCurrentRates M.! i of
        Nothing -> 1.0
        Just currentRate -> exp $ (/ 2) $ currentRate - optimalRate
      bs = mc3ReciprocalTemperatures a
      n = fromNChains $ mc3NChains $ mc3Settings a
      -- Do not change the temperature, and the prior and likelihood functions of
      -- the cold chain.
      bs' = foldl' (\xs j -> tuneBeta bs j (xi j) xs) bs [0 .. n - 2]
      coldChain = fromMHG $ V.head mhgs'
      coldPrF = priorFunction coldChain
      coldLhF = likelihoodFunction coldChain
      mhgs'' =
        V.head mhgs'
          `V.cons` V.zipWith
            (setReciprocalTemperature coldPrF coldLhF)
            (V.convert $ U.tail bs')
            (V.tail mhgs')
  return $
    if b == NormalTuningFastProposalsOnly || b == NormalTuningAllProposals
      then a {mc3MHGChains = mhgs'', mc3ReciprocalTemperatures = bs'}
      else a {mc3MHGChains = mhgs''}

mc3ResetAcceptance :: (ToJSON a) => ResetAcceptance -> MC3 a -> MC3 a
mc3ResetAcceptance x a = a'
  where
    -- 1. Reset acceptance of all chains.
    mhgs' = V.map (aResetAcceptance x) (mc3MHGChains a)
    -- 2. Reset acceptance of swaps.
    ac' = resetA x $ mc3SwapAcceptances a
    --
    a' = a {mc3MHGChains = mhgs', mc3SwapAcceptances = ac'}

mc3CleanAfterBurnIn :: (ToJSON a) => TraceLength -> MC3 a -> IO (MC3 a)
mc3CleanAfterBurnIn tl a = do
  cs' <- V.mapM (aCleanAfterBurnIn tl) cs
  pure $ a {mc3MHGChains = cs'}
  where
    cs = mc3MHGChains a

-- Information in cycle summary:
--
-- - The complete summary of the cycle of the cold chain.
--
-- - The combined acceptance rate of proposals within the hot chains.
--
-- - The temperatures of the chains and the acceptance rates of the state swaps.
mc3SummarizeCycle :: (ToJSON a) => IterationMode -> MC3 a -> BL.ByteString
mc3SummarizeCycle m a =
  BL.intercalate "\n" $
    [ "MC3: Cycle of cold chain.",
      coldMHGCycleSummary
    ]
      ++ case mAr of
        Nothing -> []
        Just ar ->
          [ "MC3: Average acceptance rate across all chains: "
              <> BB.toLazyByteString (BB.formatDouble (BB.standard 2) ar)
              <> "."
          ]
      ++ [ "MC3: Reciprocal temperatures of the chains: " <> BL.intercalate ", " bsB <> ".",
           "MC3: Summary of state swaps.",
           "MC3: The swap period is " <> swapPeriodB <> ".",
           "MC3: The state swaps are executed in random order.",
           proposalHeader,
           proposalHLine
         ]
      ++ [ summarizeProposal
             (PName $ show i ++ " <-> " ++ show (i + 1))
             (PDescription "Swap states between chains")
             (pWeight 1)
             (Just $ bs U.! (i + 1))
             PDimensionUnknown
             (acceptanceRate i swapAcceptance)
           | i <- [0 .. n - 2]
         ]
      ++ [proposalHLine]
  where
    mhgs = mc3MHGChains a
    coldMHGCycleSummary = aSummarizeCycle m $ V.head mhgs
    cs = V.map fromMHG mhgs
    -- Acceptance rates may be 'Nothing' when no proposals have been undertaken.
    -- The 'sequence' operations pull the 'Nothing's out of the inner
    -- structures.
    as = sequence $ V.map (sequence . acceptanceRates . acceptances) cs
    mVecAr = V.map (\mp -> sum mp / fromIntegral (length mp)) <$> as
    mAr = (\vec -> V.sum vec / fromIntegral (V.length vec)) <$> mVecAr
    bs = mc3ReciprocalTemperatures a
    bsB = map (BB.toLazyByteString . BB.formatDouble (BB.standard 2)) $ U.toList bs
    swapPeriod = fromSwapPeriod $ mc3SwapPeriod $ mc3Settings a
    swapPeriodB = BB.toLazyByteString $ BB.intDec swapPeriod
    swapAcceptance = mc3SwapAcceptances a
    n = fromNChains $ mc3NChains $ mc3Settings a
    proposalHLine = BL.replicate (BL.length proposalHeader) '-'

-- No extra monitors are opened.
mc3OpenMonitors :: (ToJSON a) => AnalysisName -> ExecutionMode -> MC3 a -> IO (MC3 a)
mc3OpenMonitors nm em a = do
  mhgs' <- V.imapM mhgOpenMonitors (mc3MHGChains a)
  return $ a {mc3MHGChains = mhgs'}
  where
    mhgOpenMonitors i (MHG c) = do
      m' <- mOpen pre suf em $ monitor c
      pure $ MHG c {monitor = m'}
      where
        pre = fromAnalysisName nm
        suf = printf "%02d" i

mc3ExecuteMonitors ::
  (ToJSON a) =>
  Verbosity ->
  -- Starting time.
  UTCTime ->
  -- Total number of iterations.
  Int ->
  MC3 a ->
  IO (Maybe BL.ByteString)
mc3ExecuteMonitors vb t0 iTotal a = V.head <$> V.imapM f (mc3MHGChains a)
  where
    -- The first chain honors verbosity.
    f 0 = aExecuteMonitors vb t0 iTotal
    -- All other chains are to be quiet.
    f _ = aExecuteMonitors Quiet t0 iTotal

mc3StdMonitorHeader :: (ToJSON a) => MC3 a -> BL.ByteString
mc3StdMonitorHeader = aStdMonitorHeader . V.head . mc3MHGChains

mc3CloseMonitors :: (ToJSON a) => MC3 a -> IO (MC3 a)
mc3CloseMonitors a = do
  mhgs' <- V.mapM aCloseMonitors $ mc3MHGChains a
  return $ a {mc3MHGChains = mhgs'}
