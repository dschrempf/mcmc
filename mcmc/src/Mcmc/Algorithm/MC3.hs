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
  )
where

import Data.Aeson
import qualified Data.Vector as V
import Mcmc.Algorithm
import Mcmc.Algorithm.Metropolis
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Internal.Random
import Mcmc.Monitor
import Mcmc.Proposal
import Numeric.Log
import System.Random.MWC

-- | The hotter the chain, the flatter the prior and the likelihood.
data HeatedChain a = HeatedChain
  { heatedChain :: MHG a,
    -- | Reciprocal temperature.
    _heatedBeta :: Double,
    coldPriorFunction :: PriorFunction a,
    coldLikelihoodFunction :: LikelihoodFunction a
  }

-- | A vector of heated chains.
type HeatedChains a = V.Vector (HeatedChain a)

initializeHeatedChain :: MHG a -> HeatedChain a
initializeHeatedChain m = HeatedChain m 1.0 pr lh
  where
    pr = priorFunction $ fromMHG m
    lh = likelihoodFunction $ fromMHG m

-- Be careful! When changing the temperature of the chain, the prior and
-- likelihood values of the trace are not updated! The prior and likelihood
-- values of the last link are updated.
setReciprocalTemperature :: Double -> HeatedChain a -> HeatedChain a
setReciprocalTemperature b hc = hc {heatedChain = MHG c'}
  where
    c = fromMHG $ heatedChain hc
    b' = Exp $ log b
    -- XXX: Like this, we need twice the amount of computations compared to
    -- using (pr x * lh x) ** b'. However, I don't think this is a serious
    -- problem.
    pr' = (** b') . coldPriorFunction hc
    lh' = (** b') . coldLikelihoodFunction hc
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

hcIterate :: ToJSON a => HeatedChain a -> IO (HeatedChain a)
hcIterate c = do
  a' <- aIterate a
  return $ c {heatedChain = a'}
  where
    a = heatedChain c

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
        as = V.map initializeHeatedChain cs
    return $ MC3 s (V.zipWith setReciprocalTemperature bs as) 0 0 0 g
  where
    n = mc3NChains s
    p = mc3SwapPeriod s
    -- XXX: Think about initial choice of reciprocal temperatures.
    bs = V.fromList $ map recip [1.0, 1.1 ..]

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

mc3Iterate :: ToJSON a => MC3 a -> IO (MC3 a)
mc3Iterate a = do
  -- 1. Iterate all chains
  cs' <- V.mapM hcIterate (mc3HeatedChains a)
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
