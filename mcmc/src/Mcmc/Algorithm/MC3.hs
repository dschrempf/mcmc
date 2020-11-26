{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

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
  ( HeatedChain,
    HeatedChains,
    MC3SwapType (..),
    MC3Settings (..),
    MC3 (..),
    mc3,
    mc3Save,
    mc3Load,
  )
where

import Data.Aeson
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Vector as V
import Mcmc.Algorithm
import Mcmc.Algorithm.Metropolis
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Environment
import Mcmc.Internal.Random
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log
import System.Random.MWC
import Text.Printf

-- | The hotter the chain, the flatter the prior and the likelihood.
data HeatedChain a = HeatedChain
  { heatedChain :: MHG a,
    -- | Reciprocal temperature.
    heatedBeta :: Double,
    coldPriorFunction :: PriorFunction a,
    coldLikelihoodFunction :: LikelihoodFunction a
  }

hcInitialize :: MHG a -> HeatedChain a
hcInitialize m = HeatedChain m 1.0 pr lh
  where
    pr = priorFunction $ fromMHG m
    lh = likelihoodFunction $ fromMHG m

-- TODO.
hcSave ::
  ToJSON a =>
  -- Maximum length of trace.
  Int ->
  -- Analysis name.
  String ->
  HeatedChain a ->
  IO ()
hcSave = undefined

-- TODO.
hcLoad ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  -- Analysis name.
  String ->
  IO (HeatedChain a)
hcLoad = undefined

-- Be careful! When changing the temperature of the chain, the prior and
-- likelihood values of the trace are not updated! The prior and likelihood
-- values of the last link are updated.
hcModifyReciprocalTemperature :: (Double -> Double) -> HeatedChain a -> HeatedChain a
hcModifyReciprocalTemperature f hc = hc {heatedChain = MHG c', heatedBeta = b'}
  where
    c = fromMHG $ heatedChain hc
    b = heatedBeta hc
    b' = f b
    b'LogDomain = Exp $ log b'
    -- XXX: Like this, we need twice the amount of computations compared to
    -- using (pr x * lh x) ** b'. However, I don't think this is a serious
    -- problem.
    --
    -- Also, to minimize computations, avoid modification of the reciprocal
    -- temperature for the cold chain.
    pr' = (** b'LogDomain) . coldPriorFunction hc
    lh' = (** b'LogDomain) . coldLikelihoodFunction hc
    x = state $ link c
    c' =
      c
        { priorFunction = pr',
          likelihoodFunction = lh',
          link = Link x (pr' x) (lh' x)
        }

-- Set a new link in a heated chain.
hcSetLink ::
  -- New link.
  Link a ->
  -- Chain with old link.
  HeatedChain a ->
  -- Chain with new link.
  HeatedChain a
hcSetLink l hc = hc {heatedChain = c'}
  where
    c = fromMHG $ heatedChain hc
    c' = MHG $ c {link = l}

-- Get the link from a heated chain.
hcGetLink :: HeatedChain a -> Link a
hcGetLink = link . fromMHG . heatedChain

hcGetPriorFunction :: HeatedChain a -> PriorFunction a
hcGetPriorFunction = priorFunction . fromMHG . heatedChain

hcGetLikelihoodFunction :: HeatedChain a -> LikelihoodFunction a
hcGetLikelihoodFunction = likelihoodFunction . fromMHG . heatedChain

hcApplyM :: (MHG a -> IO (MHG a)) -> HeatedChain a -> IO (HeatedChain a)
hcApplyM m c = do
  a' <- m $ heatedChain c
  return $ c {heatedChain = a'}

hcApply :: (MHG a -> MHG a) -> HeatedChain a -> HeatedChain a
hcApply f c = c {heatedChain = c'}
  where
    c' = f $ heatedChain c

-- | A vector of heated chains.
type HeatedChains a = V.Vector (HeatedChain a)

-- | Swap states between neighboring chains, or random chains.
data MC3SwapType = MC3Neighbor | MC3Random
  deriving (Eq, Read, Show)

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
    mc3HeatedChains :: HeatedChains a,
    mc3Iteration :: Int,
    -- | Number of accepted swaps.
    mc3SwapAccepted :: Int,
    -- | Number of rejected swaps.
    mc3SwapRejected :: Int,
    mc3Generator :: GenIO
  }

instance ToJSON a => Algorithm MC3 a where
  aName = const "Metropolis-coupled Markov chain Monte Carlo"
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
    let cs = V.map (mhg pr lh cc mn i0) gs
        as = V.map hcInitialize cs
        -- Do not change the prior and likelihood functions of the first chain.
        hcs = V.head as `V.cons` V.zipWith hcModifyReciprocalTemperature bs (V.tail as)
    return $ MC3 s hcs 0 0 0 g
  where
    n = mc3NChains s
    p = mc3SwapPeriod s
    -- XXX: Maybe improve initial choice of reciprocal temperatures.
    bs = V.fromList $ map (const . recip) [1.1 ..]

-- TODO.

-- | Save an MC3 algorithm.
mc3Save ::
  -- | Maximum length of trace.
  Int ->
  -- | Analysis name.
  String ->
  MC3 a ->
  IO ()
mc3Save = undefined

-- TODO.

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
mc3Load pr lh cc mn nm = undefined

-- I call the chains left and right, because it is easy to think about them as
-- being left and right. Of course, the left chain may also have a larger index
-- than the right chain.
swapWith ::
  -- Index i>=0 of left chain.
  Int ->
  -- Index j>=0, j/=i of right chain.
  Int ->
  HeatedChains a ->
  (HeatedChains a, Log Double)
swapWith i j xs
  | i < 0 = error "swapWith: Left index is negative."
  | j < 0 = error "swapWith: Right index is negative."
  | i == j = error "swapWith: Indices are equal."
  | otherwise = (xs', q)
  where
    -- Gather information from current chains.
    cl = xs V.! i
    cr = xs V.! j
    ll = hcGetLink cl
    lr = hcGetLink cr
    prl = prior ll
    prr = prior lr
    lhl = likelihood ll
    lhr = likelihood lr
    -- Swap the states.
    xl' = state lr
    xr' = state ll
    -- Compute new priors and likelihoods.
    prl' = hcGetPriorFunction cl xl'
    prr' = hcGetPriorFunction cr xr'
    lhl' = hcGetLikelihoodFunction cl xl'
    lhr' = hcGetLikelihoodFunction cr xr'
    ll' = Link xl' prl' lhl'
    lr' = Link xr' prr' lhr'
    -- Set the new links and the proposed state.
    cl' = hcSetLink ll' cl
    cr' = hcSetLink lr' cr
    xs' = xs V.// [(i, cl'), (j, cr')]
    -- Compute the Metropolis ratio.
    q = prl' * prr' * lhl' * lhr' / prl / prr / lhl / lhr

-- i is the index of left chain (0 ..).
-- j>i is the index of the right chain.
swap :: MC3SwapType -> HeatedChains a -> GenIO -> IO (HeatedChains a, Log Double)
swap st xs g =
  do
    i <- uniformR (0, n - 1) g
    j <- case st of
      MC3Neighbor -> return $ i + 1
      MC3Random -> loop i g
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
  (!y, !r) <- swap (mc3SwapType s) (mc3HeatedChains a) g
  -- 3. Accept or reject.
  accept <- mhgAccept r g
  if accept
    then do
      let !ac' = mc3SwapAccepted a + 1
      return $ a {mc3HeatedChains = y, mc3SwapAccepted = ac'}
    else do
      let !re' = mc3SwapRejected a + 1
      return $ a {mc3SwapRejected = re'}
  where
    s = mc3Settings a
    g = mc3Generator a

-- TODO: Parallel iteration.
mc3Iterate :: ToJSON a => MC3 a -> IO (MC3 a)
mc3Iterate a = do
  -- 1. Iterate all chains.
  cs' <- V.mapM (hcApplyM aIterate) (mc3HeatedChains a)
  let a' = a {mc3HeatedChains = cs'}
  -- 2. Maybe propose swap.
  let s' = mc3Settings a'
  a'' <-
    if mc3Iteration a' `mod` mc3SwapPeriod s' == 0
      then mc3ProposeSwap a'
      else return a'
  -- 3. Increment iteration counter.
  let i'' = mc3Iteration a''
  return $ a'' {mc3Iteration = succ i''}

mc3AutoTune :: ToJSON a => MC3 a -> MC3 a
mc3AutoTune a = a {mc3HeatedChains = cs''}
  where
    -- 1. Auto tune all chains.
    cs' = V.map (hcApply aAutoTune) (mc3HeatedChains a)
    -- 2. Auto tune temperatures.
    currentRate = fromIntegral (mc3SwapAccepted a) / fromIntegral (mc3SwapRejected a)
    optimalRate = getOptimalRate PDimensionUnknown
    dt = exp (currentRate - optimalRate)
    -- Do not change the prior and likelihood functions of the cold chain.
    cs'' =
      V.head cs'
        `V.cons` V.map (hcModifyReciprocalTemperature (* dt)) (V.tail cs')

mc3ResetAcceptance :: ToJSON a => MC3 a -> MC3 a
mc3ResetAcceptance a = a'
  where
    -- 1. Reset acceptance of all chains.
    cs' = V.map (hcApply aResetAcceptance) (mc3HeatedChains a)
    -- 2. Reset acceptance of swaps.
    a' = a {mc3HeatedChains = cs', mc3SwapAccepted = 0, mc3SwapRejected = 0}

-- TODO.
--
-- Information in cycle summary:
--
-- - The complete summary of the cycle of the cold chain.
--
-- - The temperatures of the chains and the acceptance rate of state swaps.
--
-- - The combined acceptance rate of proposals within the hot chains.
mc3SummarizeCycle :: MC3 a -> BL.ByteString
mc3SummarizeCycle = undefined

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
  cs' <- V.imapM f $ mc3HeatedChains a
  return $ a {mc3HeatedChains = cs'}
  where
    f i = hcApplyM $ aOpenMonitors (amendEnvironment i e)

-- TODO: Parallel execution?

-- TODO: It is a little unfortunate that we have to amend the environment every time.
mc3ExecuteMonitors :: ToJSON a => Environment -> MC3 a -> IO (Maybe BL.ByteString)
mc3ExecuteMonitors e a = V.head <$> V.imapM f (mc3HeatedChains a)
  where
    f i = aExecuteMonitors (amendEnvironment i e) . heatedChain

mc3StdMonitorHeader :: ToJSON a => MC3 a -> BL.ByteString
mc3StdMonitorHeader = aStdMonitorHeader . heatedChain . V.head . mc3HeatedChains

mc3CloseMonitors :: ToJSON a => MC3 a -> IO (MC3 a)
mc3CloseMonitors a = do
  cs' <- V.mapM (hcApplyM aCloseMonitors) $ mc3HeatedChains a
  return $ a {mc3HeatedChains = cs'}
