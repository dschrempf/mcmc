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
    MC3 (..),
    mc3,
  )
where

import Mcmc.Algorithm
import Mcmc.Algorithm.Metropolis
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Monitor
import Mcmc.Proposal
import Numeric.Log
import System.Random.MWC

-- | The hotter the chain, the flatter the prior and the likelihood.
data HeatedChain a = HeatedChain
  { heatedChain :: MHG a,
    -- | Reciprocal temperature.
    heatedBeta :: Double,
    coldPriorFunction :: PriorFunction a,
    coldLikelihoodFunction :: LikelihoodFunction a
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
setLink ::
  -- New link.
  Link a ->
  -- Chain with old link.
  HeatedChain a ->
  -- Chain with new link.
  HeatedChain a
setLink l hc = hc {heatedChain = c'}
  where
    c = fromMHG $ heatedChain hc
    c' = MHG $ c {link = l}

-- Get the link from a heated chain.
getLink :: HeatedChain a -> Link a
getLink = link . fromMHG . heatedChain

-- | The MC3 algorithm.
--
-- Also known as parallel tempering.
data MC3 a = MC3
  { -- | The first chain is the cold chain with temperature 1.0.
    mc3HeatedChains :: [HeatedChain a],
    mc3SwapPeriod :: Int
  }

-- TODO.
instance Algorithm MC3 a

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
  | otherwise = MC3 (zipWith setReciprocalTemperature bs as) p
  where
    a = mhg pr lh cc mn i0 g
    -- TODO: Think about initial choice of reciprocal temperatures.
    bs = [1.0, 1.1 ..]
    as = replicate n $ initializeHeatedChain a

-- TODO: Acceptance ratio should be 0.234, since this is a high dimensional
-- proposal. I think it makes sense to implement a dynamic acceptance ratio for
-- all proposals at this point.

-- TODO: Start implementing swap of beighboring chains. Determine the index
-- randomly.
mc3ProposeSwap :: MC3 a -> IO (MHG a)
mc3ProposeSwap = undefined
