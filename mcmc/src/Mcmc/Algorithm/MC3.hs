{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Algorithm.MC3
-- Description :  Metropolis-coupled Markov chain Monte Carlo algorithm
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 23 15:20:33 2020.
module Mcmc.Algorithm.MC3
  ( NChains (..),
    SwapPeriod (..),
    NSwaps (..),
    MC3Settings (..),
    MHGChains,
    ReciprocalTemperatures,
    MC3 (..),
    mc3,
    mc3Save,
    mc3Load,
  )
where

import Codec.Compression.GZip
import Control.Monad
import qualified Control.Monad.Parallel as P
import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Double.Conversion.ByteString as BC
import Data.List
import qualified Data.Map.Strict as M
import Data.Time
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Word
-- import Debug.Trace hiding (trace)
import Mcmc.Algorithm
import Mcmc.Algorithm.Metropolis
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Save
import Mcmc.Chain.Trace
import Mcmc.Internal.Random
import Mcmc.Internal.Shuffle
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log hiding (sum)
import System.Random.MWC

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

-- XXX: An unboxed vector could be used for 'ReciprocalTemperatures', but it is
-- complicated. Since we zip with the boxed vector 'Chains'.

-- | Vector of reciprocal temperatures.
type ReciprocalTemperatures = U.Vector Double

data SavedMC3 a = SavedMC3
  { savedMC3Settings :: MC3Settings,
    savedMC3Chains :: V.Vector (SavedChain a),
    savedMC3ReciprocalTemperatures :: ReciprocalTemperatures,
    savedMC3Iteration :: Int,
    savedMC3SwapAcceptance :: Acceptance Int,
    savedMC3Generator :: U.Vector Word32
  }
  deriving (Eq, Read, Show)

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
  mhgs <- V.mapM (fmap MHG . fromSavedChain pr lh cc mn) scs
  g <- loadGen g'
  return $ MC3 s mhgs bs i ac g

-- | The MC3 algorithm.
--
-- Also known as parallel tempering.
--
-- Geyer, C. J., Markov chain monte carlo maximum likelihood, Computing Science
-- and Statistics, Proceedings of the 23rd Symposium on the Interface, (1991).
--
-- Altekar, G., Dwarkadas, S., Huelsenbeck, J. P., & Ronquist, F., Parallel
-- metropolis coupled markov chain monte carlo for bayesian phylogenetic
-- inference, Bioinformatics, 20(3), 407–415 (2004).
data MC3 a = MC3
  { mc3Settings :: MC3Settings,
    -- | The first chain is the cold chain with temperature 1.0.
    mc3MHGChains :: MHGChains a,
    -- | Vector of reciprocal temperatures.
    mc3ReciprocalTemperatures :: ReciprocalTemperatures,
    mc3Iteration :: Int,
    -- | Number of accepted and rejected swaps.
    mc3SwapAcceptance :: Acceptance Int,
    mc3Generator :: GenIO
  }

instance ToJSON a => Algorithm (MC3 a) where
  aName = const "Metropolis-coupled Markov chain Monte Carlo (MC3)"
  aIteration = mc3Iteration
  aIterate = mc3Iterate
  aAutoTune = mc3AutoTune
  aResetAcceptance = mc3ResetAcceptance
  aSummarizeCycle = mc3SummarizeCycle
  aOpenMonitors = mc3OpenMonitors
  aExecuteMonitors = mc3ExecuteMonitors
  aStdMonitorHeader = mc3StdMonitorHeader
  aCloseMonitors = mc3CloseMonitors
  aSave = mc3Save

--  The prior and likelihood values of the current link are updated.
--
-- Be careful. The trace is not changed! In particular, the prior and likelihood
-- values are not updated for any link of the trace.
setReciprocalTemperature ::
  -- Cold prior function.
  PriorFunction a ->
  -- Cold likelihood function.
  LikelihoodFunction a ->
  -- New reciprocal temperature.
  Double ->
  MHG a ->
  MHG a
setReciprocalTemperature prf lhf beta a =
  MHG $
    c
      { priorFunction = prf',
        likelihoodFunction = lhf',
        link = Link x (prf' x) (lhf' x)
      }
  where
    c = fromMHG a
    b' = Exp $ log beta
    -- Like this, we need twice the amount of computations compared to taking
    -- the power after calculating the posterior (pr x * lh x) ** b'. However, I
    -- don't think this is a serious problem.
    --
    -- Also, to minimize computations, avoid modification of the reciprocal
    -- temperature for the cold chain.
    prf' = (** b') . prf
    lhf' = (** b') . lhf
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
  -- Do not temper with the first chain.
  | i == 0 = return a
  | otherwise = do
    -- We have to reset the trace, since it is not set by 'setReciprocalTemperature'.
      t' <- pushT l t
      return $ MHG $ c {chainId = i + 1, trace = t'}
  where
    a' = setReciprocalTemperature prf lhf beta a
    c = fromMHG a'
    l = link c
    t = trace c

-- TODO: Splitmix. Initialization of the MC3 algorithm is an IO action because
-- the generators have to be split. And also because of the mutable trace.

-- | Initialize an MC3 algorithm with a given number of chains.
--
-- Call 'error' if:
--
-- - The number of chains is one or lower.
--
-- - The swap period is zero or negative.
mc3 ::
  MC3Settings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  a ->
  GenIO ->
  IO (MC3 a)
mc3 s pr lh cc mn i0 g
  | n < 2 = error "mc3: The number of chains must be two or larger."
  | sp < 1 = error "mc3: The swap period must be strictly positive."
  | sn < 1 || sn > n - 1 = error "mc3: The number of swaps must be in [1, NChains - 1]."
  | otherwise = do
    -- Split random number generators.
    gs <- V.fromList <$> splitGen n g
    cs <- V.mapM (mhg pr lh cc mn i0) gs
    hcs <- V.izipWithM (initMHG pr lh) (V.convert bs) cs
    return $ MC3 s hcs bs 0 (emptyA [0 .. n - 2]) g
  where
    n = fromNChains $ mc3NChains s
    sp = fromSwapPeriod $ mc3SwapPeriod s
    sn = fromNSwaps $ mc3NSwaps s
    -- XXX: Maybe improve initial choice of reciprocal temperatures.
    --
    -- Have to 'take n' elements, because vectors are not as lazy as lists.
    bs = U.fromList $ take n $ iterate (* 0.9) 1.0

mc3Fn :: AnalysisName -> FilePath
mc3Fn (AnalysisName nm) = nm ++ ".mc3"

-- | Save an MC3 algorithm.
mc3Save ::
  ToJSON a =>
  AnalysisName ->
  MC3 a ->
  IO ()
mc3Save nm a = do
  savedMC3 <- toSavedMC3 a
  BL.writeFile (mc3Fn nm) $ compress $ encode savedMC3

-- | Load an MC3 algorithm.
mc3Load ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  AnalysisName ->
  IO (MC3 a)
mc3Load pr lh cc mn nm = do
  savedMC3 <- eitherDecode . decompress <$> BL.readFile (mc3Fn nm)
  either error (fromSavedMC3 pr lh cc mn) savedMC3

-- I call the chains left and right, because it is easy to think about them as
-- being left and right. Of course, the left chain may also have a larger index
-- than the right chain.
swapWith ::
  -- Index i>=0 of left chain.
  Int ->
  -- Index j>=0, j/=i of right chain.
  Int ->
  MHGChains a ->
  (MHGChains a, Log Double)
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
    ll' = Link xl' prl' lhl'
    lr' = Link xr' prr' lhr'
    -- Set the new links and the proposed state.
    cl' = cl {link = ll'}
    cr' = cr {link = lr'}
    xs' = xs V.// [(i, MHG cl'), (j, MHG cr')]
    -- Compute the Metropolis ratio.
    q = prl' * prr' * lhl' * lhr' / prl / prr / lhl / lhr

mc3ProposeSwap :: MC3 a -> Int -> IO (MC3 a)
mc3ProposeSwap a i = do
  -- 1. Sample new state and get the Metropolis ratio.
  let (!y, !r) = swapWith i (i + 1) $ mc3MHGChains a
  -- 2. Accept or reject.
  accept <- mhgAccept r g
  if accept
    then do
      let !ac' = pushA i True (mc3SwapAcceptance a)
      return $ a {mc3MHGChains = y, mc3SwapAcceptance = ac'}
    else do
      let !ac' = pushA i False (mc3SwapAcceptance a)
      return $ a {mc3SwapAcceptance = ac'}
  where
    g = mc3Generator a

-- TODO: Splimix. 'mc3Iterate' is actually not parallel, but concurrent because
-- of the IO constraint. Use pure parallel code when we have a pure generator.
mc3Iterate ::
  ToJSON a =>
  ParallelizationMode ->
  MC3 a ->
  IO (MC3 a)
mc3Iterate pm a = do
  -- 1. Maybe propose swap. A swap has to be proposed first, because the traces
  -- are automatically updated at step 2.
  let s = mc3Settings a
  a' <-
    if mc3Iteration a `mod` fromSwapPeriod (mc3SwapPeriod s) == 0
      then do
        let n = V.length $ mc3MHGChains a
            is = [0 .. n - 2]
            ns = fromNSwaps $ mc3NSwaps $ mc3Settings a
        is' <- shuffle is $ mc3Generator a
        foldM mc3ProposeSwap a (take ns is')
      else return a
  -- 2. Iterate all chains and increment iteration.
  --
  -- In any case, the MHG algorithm only gets one capability.
  mhgs <- case pm of
    Sequential -> V.mapM (aIterate pm) (mc3MHGChains a')
    Parallel ->
      -- Use 'forkIO'.
      V.fromList <$> P.mapM (aIterate pm) (V.toList (mc3MHGChains a'))
  let i = mc3Iteration a'
  return $ a' {mc3MHGChains = mhgs, mc3Iteration = succ i}

tuneBeta ::
  -- The old reciprocal temperatures are needed to retrieve the old ratios.
  ReciprocalTemperatures ->
  -- Index i of left chain. Change the reciprocal temperature of chain (i+1).
  Int ->
  -- Exponent xi of the reciprocal temperature ratio.
  Double ->
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

mc3AutoTune :: ToJSON a => MC3 a -> MC3 a
mc3AutoTune a = a {mc3MHGChains = mhgs'', mc3ReciprocalTemperatures = bs'}
  where
    mhgs = mc3MHGChains a
    -- 1. Auto tune all chains.
    mhgs' = V.map aAutoTune mhgs
    -- 2. Auto tune temperatures.
    coldChain = fromMHG $ V.head mhgs
    coldPrF = priorFunction coldChain
    coldLhF = likelihoodFunction coldChain
    optimalRate = getOptimalRate PDimensionUnknown
    currentRates = acceptanceRates $ mc3SwapAcceptance a
    -- We assume that the acceptance rate of state swaps between two chains is
    -- roughly proportional to the ratio of the temperatures of the chains.
    -- Hence, we focus on reciprocal temperature ratios, which is the same.
    -- Also, by working with ratios in (0,1) of neighboring chains, we ensure
    -- the monotonicity of the reciprocal temperatures.
    --
    -- The factor (1/2) was determined by a few tests and is otherwise
    -- absolutely arbitrary.
    xi i = exp $ (/ 2) $ (currentRates M.! i) - optimalRate
    bs = mc3ReciprocalTemperatures a
    n = fromNChains $ mc3NChains $ mc3Settings a
    -- Do not change the temperature, and the prior and likelihood functions of
    -- the cold chain.
    bs' = foldl' (\xs j -> tuneBeta bs j (xi j) xs) bs [0 .. n - 2]
    mhgs'' =
      V.head mhgs'
        `V.cons` V.zipWith
          (setReciprocalTemperature coldPrF coldLhF)
          (V.convert $ U.tail bs')
          (V.tail mhgs')

mc3ResetAcceptance :: ToJSON a => MC3 a -> MC3 a
mc3ResetAcceptance a = a'
  where
    -- 1. Reset acceptance of all chains.
    mhgs' = V.map aResetAcceptance (mc3MHGChains a)
    -- 2. Reset acceptance of swaps.
    ac = mc3SwapAcceptance a
    a' = a {mc3MHGChains = mhgs', mc3SwapAcceptance = resetA ac}

-- Information in cycle summary:
--
-- - The complete summary of the cycle of the cold chain.
--
-- - The combined acceptance rate of proposals within the hot chains.
--
-- - The temperatures of the chains and the acceptance rates of the state swaps.
mc3SummarizeCycle :: ToJSON a => MC3 a -> BL.ByteString
mc3SummarizeCycle a =
  BL.intercalate "\n" $
    [ "MC3: Cycle of cold chain.",
      coldMHGCycleSummary
    ]
      ++ [ "MC3: Average acceptance rate across all chains: " <> BL.fromStrict (BC.toFixed 2 ar)
           | not $ isNaN ar
         ]
      ++ [ "MC3: Reciprocal temperatures of the chains: " <> BL.intercalate ", " bsB <> ".",
           "MC3: Summary of state swaps. The swap period is " <> swapPeriodB <> ".",
           proposalHeader,
           proposalHLine
         ]
      ++ [ summarizeProposal
             (PName $ show i ++ " <-> " ++ show (i + 1))
             (PDescription "Swap states between chains")
             (PWeight 1)
             (Just $ bs U.! (i + 1))
             PDimensionUnknown
             (acceptanceRate i swapAcceptance)
           | i <- [0 .. n - 2]
         ]
      ++ [proposalHLine]
  where
    mhgs = mc3MHGChains a
    coldMHGCycleSummary = aSummarizeCycle $ V.head mhgs
    cs = V.map fromMHG mhgs
    as = V.map (acceptanceRates . acceptance) cs
    vAr = V.map (\m -> sum m / fromIntegral (length m)) as
    ar = V.sum vAr / fromIntegral (V.length vAr)
    bs = mc3ReciprocalTemperatures a
    bsB = map (BL.fromStrict . BC.toFixed 2) $ U.toList bs
    swapPeriod = fromSwapPeriod $ mc3SwapPeriod $ mc3Settings a
    swapPeriodB = BB.toLazyByteString $ BB.intDec swapPeriod
    swapAcceptance = mc3SwapAcceptance a
    n = fromNChains $ mc3NChains $ mc3Settings a

-- See 'amendEnvironment'.
--
-- No extra monitors are opened.
mc3OpenMonitors :: ToJSON a => AnalysisName -> ExecutionMode -> MC3 a -> IO (MC3 a)
mc3OpenMonitors nm em a = do
  mhgs' <- V.mapM (aOpenMonitors nm em) (mc3MHGChains a)
  return $ a {mc3MHGChains = mhgs'}

mc3ExecuteMonitors ::
  ToJSON a =>
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

mc3StdMonitorHeader :: ToJSON a => MC3 a -> BL.ByteString
mc3StdMonitorHeader = aStdMonitorHeader . V.head . mc3MHGChains

mc3CloseMonitors :: ToJSON a => MC3 a -> IO (MC3 a)
mc3CloseMonitors a = do
  mhgs' <- V.mapM aCloseMonitors $ mc3MHGChains a
  return $ a {mc3MHGChains = mhgs'}
