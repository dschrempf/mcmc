-- |
-- Module      :  Mcmc.Algorithm.Mc3
-- Description :  Metropolis-coupled Markov chain Monte Carlo algorithm
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 23 15:20:33 2020.
module Mcmc.Algorithm.Mc3
  ( HeatedChain,
    MC3 (..),
    mc3,
  )
where

import Mcmc.Algorithm.Metropolis
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Monitor
import Mcmc.Proposal
import Numeric.Log
import System.Random.MWC

-- | The hotter the chain, the flatter the prior and the likelihood.
data HeatedChain a = HeatedChain
  { _heatedChain :: MHG a,
    -- | Reciprocal temperature.
    _heatedBeta :: Double,
    _coldPriorFunction :: PriorFunction a,
    _coldLikelihoodFunction :: LikelihoodFunction a
  }

initializeHeatedChain :: MHG a -> HeatedChain a
initializeHeatedChain m = HeatedChain m 1.0 pr lh
  where
    pr = priorFunction $ fromMHG m
    lh = likelihoodFunction $ fromMHG m

-- Be careful! When changing the temperature of the chain, the prior and
-- likelihood values of the trace are not updated! The prior and likelihood
-- values of the last link are updated.
setReciprocalTemperature :: Double -> HeatedChain a -> HeatedChain a
setReciprocalTemperature b hc = hc {_heatedChain = MHG c'}
  where
    c = fromMHG $ _heatedChain hc
    b' = Exp $ log b
    pr' = (** b') . _coldPriorFunction hc
    lh' = (** b') . _coldLikelihoodFunction hc
    x = state $ link c
    c' =
      c
        { priorFunction = pr',
          likelihoodFunction = lh',
          link = Link x (pr' x) (lh' x)
        }

-- | The MC3 algorithm.
--
-- Also known as parallel tempering.
data MC3 a = MC3
  { mc3ColdChain :: MHG a,
    mc3HeatedChains :: [HeatedChain a],
    mc3SwapPeriod :: Int
  }

-- | Initialize an MC3 algorithm with a given number of chains.
mc3 ::
  -- | The number of chains.
  Int ->
  -- | Period to propose temperature swaps.
  Int ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  a ->
  GenIO ->
  MC3 a
mc3 n p pr lh cc mn i0 g
  | n < 2 = error "mc3: The number of chains must be two or larger."
  | p < 1 = error "mc3: The swap period must be strictly positive."
  | otherwise = MC3 a (zipWith setReciprocalTemperature bs as) p
  where
    a = mhg pr lh cc mn i0 g
    -- TODO: Think about initial choice of reciprocal temperatures.
    bs = [1.1, 1.2 ..]
    as = replicate (n - 1) $ initializeHeatedChain a
