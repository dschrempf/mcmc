-- |
-- Module      :  Mcmc.Tree.Prior.RelaxedClock
-- Description :  Relaxed clock models
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Sep 10 13:53:10 2020.
module Mcmc.Tree.Prior.RelaxedClock
  (
    uncorrelatedGamma,
    uncorrelatedGammaNoStem,
    whiteNoise,
    whiteNoiseNoStem,
  )
where

import ELynx.Tree
import Mcmc.Prior
import Mcmc.Tree.Prior.Branch
import Numeric.Log

-- | Uncorrelated gamma model.
--
-- The rates are distributed according to a gamma distribution with given shape
-- and scale.
--
-- For a version that ignores the root branch, see 'uncorrelatedGamma''.
uncorrelatedGamma :: Double -> Double -> Tree Double a -> Log Double
uncorrelatedGamma k th = branchesWith (gamma k th)

-- | See 'uncorrelatedGamma' but ignore the stem.
uncorrelatedGammaNoStem :: Double -> Double -> Tree Double a -> Log Double
uncorrelatedGammaNoStem k th = branchesWith' (gamma k th)

-- | White noise model.
--
-- The rates are distributed according to a white noise process with given
-- variance.
--
-- The time tree normalized to height 1.0 has to be given because long branches
-- are expected to have a distribution of rates with a lower variance than short
-- branches.
--
-- For a version that ignores the root branch, see 'whiteNoise''.
--
-- Gives unexpected results if the topologies do not match.
whiteNoise :: Double -> Tree Double a -> Tree Double a -> Log Double
whiteNoise v t r = gamma k (1 / k) (branch r) * whiteNoiseNoStem v t r
  where
    k = branch t / v

-- | See 'whiteNoise' but ignore the stem.
whiteNoiseNoStem :: Double -> Tree Double a -> Tree Double a -> Log Double
whiteNoiseNoStem v (Node _ _ ts) (Node _ _ rs) = product $ zipWith (whiteNoise v) ts rs
