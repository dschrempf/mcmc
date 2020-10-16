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
  ( gammaDirichlet,
    uncorrelatedGamma,
    uncorrelatedGammaNoStem,
    whiteNoise,
    whiteNoiseNoStem,
  )
where

import ELynx.Tree
import Mcmc.Prior
import Mcmc.Tree.Prior.Branch
import Numeric.Log hiding (sum)
import Numeric.SpecFunctions

-- Tolerance.
eps :: Double
eps = 1e-12

dirichletSymmetric :: Double -> [Double] -> Log Double
dirichletSymmetric alpha xs =
  if abs (sum xs - 1.0) > eps
  then 0
  else Exp $ logDenominator - logNominator + sum xsPow
  where
    n = length xs
    logNominator = fromIntegral n * logGamma alpha
    logDenominator = logGamma (fromIntegral n * alpha)
    xsPow = map (\x -> log $ x ** (alpha - 1.0)) xs

-- | Gamma Dirichlet prior.
--
-- Use a gamma distribution with shape parameter \(\alpha_{\mu}\) and scale
-- parameter \(\beta_{\mu}\) as prior for the mean rate \(\bar{\mu}\) of \(L\)
-- partitions. Use a symmetric Dirichlet distribution prior with given
-- concentration parameter \(\alpha\) to distribute the total rate \(\bar{\mu}\)
-- across the \(L\) partitions.
--
-- Usage:
--
-- @
--   gammaDirichlet alphaMu betaMu alphaDirichlet muMean relativeRates
-- @
--
-- The actual rate of a partition \(i\) is then calculated as
-- \[
-- \mu_i = x_i * L * \bar{\mu},
-- \]
-- where \(x_i\) is the relative rate of partition \(i\).
--
-- See Dos Reis, M., Zhu, T., & Yang, Z., The impact of the rate prior on
-- bayesian estimation of divergence times with multiple loci, Systematic
-- Biology, 63(4), 555–565 (2014). http://dx.doi.org/10.1093/sysbio/syu020.
--
-- Note that here, the SCALE and not the RATE is used as parameter of the gamma
-- distribution (in contrast to the cited publication above).
--
-- Return a probability of zero if the relative rates do not sum to 1.0 (with
-- tolerance 1e-12).
gammaDirichlet :: Double -> Double -> Double -> Double -> [Double] -> Log Double
gammaDirichlet alphaMu betaMu alpha muMean xs = muPrior * dirichletSymmetric alpha xs
  where muPrior = gamma alphaMu betaMu muMean

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
uncorrelatedGammaNoStem k th = branchesWithNoStem (gamma k th)

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
