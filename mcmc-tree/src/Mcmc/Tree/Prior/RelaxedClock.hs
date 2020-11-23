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
    uncorrelatedLogNormal,
    whiteNoise,
  )
where

-- XXX: Auto correlated UGAM model: https://doi.org/10.1111/mec.12953.

import qualified Data.Vector.Unboxed as V
import ELynx.Tree
import Mcmc.Prior
import Mcmc.Tree.Prior.Branch
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)
import Numeric.MathFunctions.Constants
import Statistics.Distribution.Dirichlet

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
--
-- Call 'error' if:
--
-- - The \(\alpha\) parameter is negative or zero.
--
-- - The number of partitions is smaller than two.
gammaDirichlet :: Double -> Double -> Double -> Double -> V.Vector Double -> Log Double
gammaDirichlet alphaMu betaMu alpha muMean xs = muPrior * dirichletDensitySymmetric ddSym xs
  where
    muPrior = gamma alphaMu betaMu muMean
    ddSym = either error id $ dirichletDistributionSymmetric (V.length xs) alpha

-- | Uncorrelated gamma model.
--
-- The rates are distributed according to a gamma distribution with given shape
-- and scale.
uncorrelatedGamma :: HandleStem -> Double -> Double -> Tree Length a -> Log Double
uncorrelatedGamma s k th = branchesWith s (gamma k th)

-- A variant of the log normal distribution. See Yang 2006, equation (7.23).
logNormal :: Double -> Double -> Double -> Log Double
logNormal mu var r = Exp $ negate t - e
  where
    t = m_ln_sqrt_2_pi + log (r * sqrt var)
    a = recip $ 2 * var
    b = log (r / mu) + 0.5 * var
    e = a * (b ** 2)

-- | Uncorrelated log normal model.
--
-- The rates are distributed according to a log normal distribution with given
-- mean and variance.
--
-- See Computational Molecular Evolution (Yang, 2006), Section 7.4.
uncorrelatedLogNormal :: HandleStem -> Double -> Double -> Tree Length a -> Log Double
uncorrelatedLogNormal s mu var = branchesWith s (logNormal mu var)

-- | Auto-correlated white noise model.
--
-- The rates are distributed according to a white noise process with given
-- variance.
--
-- The time tree normalized to height 1.0 has to be given because long branches
-- are expected to have a distribution of rates with a lower variance than short
-- branches.
--
-- Gives unexpected results if the topologies do not match.
whiteNoise :: HandleStem -> Double -> Tree Length a -> Tree Length a -> Log Double
whiteNoise WithStem v t r = gamma k (1 / k) (fromLength $ branch r) * whiteNoise WithoutStem v t r
  where
    k = fromLength (branch t) / v
whiteNoise WithoutStem v t r = product $ zipWith (whiteNoise WithStem v) (forest t) (forest r)
