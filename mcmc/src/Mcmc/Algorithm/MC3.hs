{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
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
  ( MC3SwapType (..),
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
import qualified Control.Monad.Parallel as P
import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Double.Conversion.ByteString as BC
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
import Mcmc.Environment
import Mcmc.Internal.Random
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log hiding (sum)
import System.Random.MWC
import Text.Printf

-- | Swap states between neighboring chains, or random chains.
data MC3SwapType = MC3SwapNeighbors | MC3SwapRandom
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''MC3SwapType)

-- | MC3 settings.
data MC3Settings = MC3Settings
  { -- | The number of chains has to be larger equal two.
    mc3NChains :: Int,
    -- | The period of proposing state swaps between chains has to be strictly
    -- positive.
    mc3SwapPeriod :: Int,
    mc3SwapType :: MC3SwapType
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''MC3Settings)

-- | Vector of MHG chains.
type MHGChains a = V.Vector (MHG a)

-- XXX: An unboxed vector could be used for 'ReciprocalTemperatures', but it is
-- complicated. Since we zip with the boxed vector 'Chains'.

-- | Vector of reciprocal temperatures.
type ReciprocalTemperatures = U.Vector Double

data MC3Saved a = MC3Saved
  { mc3SavedSettings :: MC3Settings,
    mc3SavedChains :: V.Vector (SavedChain a),
    mc3SavedReciprocalTemperatures :: ReciprocalTemperatures,
    mc3SavedIteration :: Int,
    mc3SavedSwapAccepted :: Int,
    mc3SavedSwapRejected :: Int,
    mc3SavedGenerator :: U.Vector Word32
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''MC3Saved)

toMC3Saved ::
  Int ->
  MC3 a ->
  MC3Saved a
toMC3Saved n (MC3 s mhgs bs i ac re g) = MC3Saved s scs bs i ac re g'
  where
    scs = V.map (toSavedChain n . fromMHG) mhgs
    g' = saveGen g

fromMC3Saved ::
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  MC3Saved a ->
  MC3 a
fromMC3Saved pr lh cc mn (MC3Saved s scs bs i ac re g') =
  MC3 s mhgs bs i ac re g
  where
    mhgs = V.map (MHG . fromSavedChain pr lh cc mn) scs
    g = loadGen g'

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
    -- | Number of accepted swaps.
    mc3SwapAccepted :: Int,
    -- | Number of rejected swaps.
    mc3SwapRejected :: Int,
    mc3Generator :: GenIO
  }

instance ToJSON a => Algorithm MC3 a where
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

-- Be careful! When changing the temperature of the chain, the prior and
-- likelihood values of the trace are not updated! The prior and likelihood
-- values of the last link are updated.
setReciprocalTemperature ::
  -- Cold prior function.
  PriorFunction a ->
  -- Cold likelihood function.
  LikelihoodFunction a ->
  -- New reciprocal temperature.
  Double ->
  MHG a ->
  MHG a
setReciprocalTemperature prf lhf b a =
  MHG $
    c
      { priorFunction = prf',
        likelihoodFunction = lhf',
        link = Link x (prf' x) (lhf' x)
      }
  where
    c = fromMHG a
    b' = Exp $ log b
    -- XXX: Like this, we need twice the amount of computations compared to
    -- using (pr x * lh x) ** b'. However, I don't think this is a serious
    -- problem.
    --
    -- Also, to minimize computations, avoid modification of the reciprocal
    -- temperature for the cold chain.
    prf' = (** b') . prf
    lhf' = (** b') . lhf
    x = state $ link c

-- TODO: Splitmix. Initialization of the MC3 algorithm is an IO action because
-- the generators have to be split.

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
  | p < 1 = error "mc3: The swap period must be strictly positive."
  | otherwise = do
    -- Split random number generators.
    gs <- V.fromList <$> splitGen n g
    -- TODO: Check if split generators are all unique and different from g.
    -- TODO: How many chains are created?
    let cs = V.map (mhg pr lh cc mn i0) gs
        -- Do not change the prior and likelihood functions of the first chain.
        hcs =
          V.head cs
            `V.cons` V.map
              resetTrace
              (V.zipWith (setReciprocalTemperature pr lh) (V.convert $ U.tail bs) (V.tail cs))
    return $ MC3 s hcs bs 0 0 0 g
  where
    n = mc3NChains s
    p = mc3SwapPeriod s
    -- XXX: Maybe improve initial choice of reciprocal temperatures.
    --
    -- Have to 'take n' elements, because vectors are not as lazy as lists.
    bs = U.fromList $ take n $ iterate (* 0.9) 1.0
    -- XXX: We have to reset the trace, since it is not set by 'setReciprocalTemperature'.
    resetTrace a =
      let c = fromMHG a
          l = link c
       in MHG $ c {trace = singletonT l}

mc3Fn :: String -> FilePath
mc3Fn nm = nm ++ ".mc3chain"

-- | Save an MC3 algorithm.
mc3Save ::
  ToJSON a =>
  -- | Maximum length of trace.
  Int ->
  -- | Analysis name.
  String ->
  MC3 a ->
  IO ()
mc3Save n nm a = BL.writeFile (mc3Fn nm) $ compress $ encode $ toMC3Saved n a

-- | Load an MC3 algorithm.
mc3Load ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  -- | Analysis name.
  String ->
  IO (MC3 a)
mc3Load pr lh cc mn nm =
  either error (fromMC3Saved pr lh cc mn)
    . eitherDecode
    . decompress
    <$> BL.readFile (mc3Fn nm)

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

-- i is the index of left chain (0 ..).
-- j>i is the index of the right chain.
swap :: MC3SwapType -> MHGChains a -> GenIO -> IO (MHGChains a, Log Double)
swap MC3SwapNeighbors xs g = do
  i <- uniformR (0, V.length xs - 2) g
  return $ swapWith i (i + 1) xs
swap MC3SwapRandom xs g = do
  i <- uniformR (0, n - 1) g
  j <- loop i g
  return $ swapWith i j xs
  where
    n = V.length xs
    -- Sample j until it is not i. This is an infinite loop if n == 1, but we
    -- ensure n > 1 in 'mc3'.
    loop l gen = do
      m <- uniformR (0, n - 1) gen
      if l == m then loop l gen else return m

mc3ProposeSwap :: MC3 a -> IO (MC3 a)
mc3ProposeSwap a = do
  -- 1. Sample new state and get the Metropolis ratio.
  (!y, !r) <- swap (mc3SwapType s) (mc3MHGChains a) g
  -- 2. Accept or reject.
  accept <- mhgAccept r g
  if accept
    then do
      let !ac' = mc3SwapAccepted a + 1
      return $ a {mc3MHGChains = y, mc3SwapAccepted = ac'}
    else do
      let !re' = mc3SwapRejected a + 1
      return $ a {mc3SwapRejected = re'}
  where
    s = mc3Settings a
    g = mc3Generator a

-- TODO: Splimix. 'mc3Iterate' is actually not parallel, but concurrent because
-- of the IO constraint. Use pure parallel code when we have a pure generator.
--
-- Parallel execution of the chains is only beneficial when the computations are
-- expensive. For all of my tests, this is not the case. However, for
-- phylogenetic dating the speed up is substantial.
mc3Iterate ::
  ToJSON a =>
  -- Number of capabilities.
  Int ->
  MC3 a ->
  IO (MC3 a)
mc3Iterate c a
  | c > 0 = do
    -- 1. Maybe propose swap. A swap has to be propose first, because the traces
    -- are automatically updated at step 2.
    let s = mc3Settings a
    a' <-
      if mc3Iteration a `mod` mc3SwapPeriod s == 0
        then mc3ProposeSwap a
        else return a
    -- 2. Iterate all chains and increment iteration.
    --
    -- In any case, the MHG algorithm only gets one capability.
    mhgs <- if c > 1
      -- Use 'forkIO'.
      then V.fromList <$> P.mapM (aIterate 1) (V.toList (mc3MHGChains a'))
      else V.mapM (aIterate 1) (mc3MHGChains a')
    let i = mc3Iteration a'
    return $ a' {mc3MHGChains = mhgs, mc3Iteration = succ i}
  | otherwise = error " mc3Iterate: Number of capabilities is zero or negative."

mc3GetAcceptanceRate :: MC3 a -> Double
mc3GetAcceptanceRate a = swapAc / (swapAc + swapRe)
  where
    swapAc = fromIntegral (mc3SwapAccepted a)
    swapRe = fromIntegral (mc3SwapRejected a)

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
    currentRate = mc3GetAcceptanceRate a
    -- The factor (1/4) was determined by a few tests and is otherwise
    -- absolutely arbitrary.
    xi = exp $ (/ 4) $ currentRate - optimalRate
    bs = mc3ReciprocalTemperatures a
    -- The reciprocal temperatures are changed differently for each chain so
    -- that the monotonicity of the reciprocal temperatures is retained.
    bf b = b ** xi
    -- Do not change the temperature, and the prior and likelihood functions of
    -- the cold chain.
    bs' = U.head bs `U.cons` U.map bf (U.tail bs)
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
    a' = a {mc3MHGChains = mhgs', mc3SwapAccepted = 0, mc3SwapRejected = 0}

-- TODO: Improve output.
--
-- Information in cycle summary:
--
-- - The complete summary of the cycle of the cold chain.
--
-- - The temperatures of the chains and the acceptance rate of state swaps.
--
-- - The combined acceptance rate of proposals within the hot chains.
mc3SummarizeCycle :: ToJSON a => MC3 a -> BL.ByteString
mc3SummarizeCycle a =
  BL.intercalate "\n" $
    [ "MC3: Cycle of cold chain.",
      coldMHGCycleSummary,
      "MC3: Average acceptance rate across all chains: " <> BL.fromStrict (BC.toFixed 2 ar),
      "MC3: Reciprocal temperatures of the chains: " <> BL.unwords bs
    ]
      ++ [ BB.toLazyByteString $
             BB.lazyByteString "MC3: State swaps: "
               <> BB.intDec swapAc
               <> BB.lazyByteString " accepted, "
               <> BB.intDec swapTot
               <> BB.lazyByteString " total, "
               <> BB.byteString (BC.toFixed 2 swapAr)
               <> BB.lazyByteString " rate."
           | swapTot /= 0
         ]
  where
    mhgs = mc3MHGChains a
    coldMHGCycleSummary = aSummarizeCycle $ V.head mhgs
    cs = V.map fromMHG mhgs
    as = V.map (acceptanceRates . acceptance) cs
    vAr = V.map (\m -> sum m / fromIntegral (length m)) as
    ar = V.sum vAr / fromIntegral (V.length vAr)
    bs = map (BL.fromStrict . BC.toFixed 2) $ U.toList $ mc3ReciprocalTemperatures a
    swapAc = mc3SwapAccepted a
    swapTot = mc3SwapRejected a + swapAc
    swapAr = mc3GetAcceptanceRate a

-- Amend the environment of chain i:
--
-- - Add a number to the analysis name.
--
-- - Do not temper with verbosity of cold chain, but set verbosity of all other
-- - chains to 'Quiet'.
amendEnvironment :: Int -> Environment -> Environment
amendEnvironment i e = e {settings = s''}
  where
    s = settings e
    n = sAnalysisName s
    suf = printf "%02d" i
    s' = s {sAnalysisName = n ++ suf}
    s'' = if i == 0 then s' else s' {sVerbosity = Quiet}

-- See 'amendEnvironment'.
--
-- No extra monitors are opened.
mc3OpenMonitors :: ToJSON a => Environment -> MC3 a -> IO (MC3 a)
mc3OpenMonitors e a = do
  mhgs' <- V.imapM f $ mc3MHGChains a
  return $ a {mc3MHGChains = mhgs'}
  where
    f j = aOpenMonitors (amendEnvironment j e)

-- TODO: Parallel execution?

-- TODO: It is a little unfortunate that we have to amend the environment every time.
mc3ExecuteMonitors :: ToJSON a => Environment -> MC3 a -> IO (Maybe BL.ByteString)
mc3ExecuteMonitors e a = V.head <$> V.imapM f (mc3MHGChains a)
  where
    f i = aExecuteMonitors (amendEnvironment i e)

mc3StdMonitorHeader :: ToJSON a => MC3 a -> BL.ByteString
mc3StdMonitorHeader = aStdMonitorHeader . V.head . mc3MHGChains

mc3CloseMonitors :: ToJSON a => MC3 a -> IO (MC3 a)
mc3CloseMonitors a = do
  mhgs' <- V.mapM aCloseMonitors $ mc3MHGChains a
  return $ a {mc3MHGChains = mhgs'}
