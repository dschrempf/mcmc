-- |
-- Module      :  Mcmc.Tree.Prior.Branch.RelaxedClock
-- Description :  Relaxed clock models
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Sep 10 13:53:10 2020.
--
-- Comparison of molecular clock models:
--
-- - Lepage, T., Bryant, D., Philippe, H., & Lartillot, N., A general comparison
--   of relaxed molecular clock models, Molecular Biology and Evolution, 24(12),
--   2669–2680 (2007). http://dx.doi.org/10.1093/molbev/msm193.
--
-- - Ho, S. Y. W., & Duchêne, S., Molecular-clock methods for estimating
--   evolutionary rates and timescales, Molecular Ecology, 23(24), 5947–5965
--   (2014). http://dx.doi.org/10.1111/mec.12953.
module Mcmc.Tree.Prior.Branch.RelaxedClock
  ( -- * Models combining relaxed molecular clocks
    gammaDirichlet,

    -- * Uncorrelated models
    uncorrelatedGamma,
    uncorrelatedLogNormal,
    whiteNoise,

    -- * Auto-correlated models
    autocorrelatedGamma,
    autocorrelatedLogNormal,
  )
where

import Data.Maybe
import qualified Data.Vector.Unboxed as V
import ELynx.Tree
import Mcmc.Chain.Chain
import Mcmc.Prior
import Mcmc.Statistics.Types
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
gammaDirichlet :: Shape -> Scale -> Double -> Mean -> PriorFunction (V.Vector Double)
gammaDirichlet alphaMu betaMu alpha muMean xs = muPrior * dirichletDensitySymmetric ddSym xs
  where
    muPrior = gamma alphaMu betaMu muMean
    ddSym = either error id $ dirichletDistributionSymmetric (V.length xs) alpha

-- | Uncorrelated gamma model.
--
-- The rates are distributed according to a gamma distribution with given mean
-- and variance.
--
-- NOTE: For convenience, the mean and variance are used as parameters for this
-- relaxed molecular clock model. They are used to calculate the shape and the
-- scale of the underlying gamma distribution.
uncorrelatedGamma :: HandleLayer -> Mean -> Variance -> PriorFunction (Tree Length a)
uncorrelatedGamma hl m v = branchesWith hl (gamma k th . fromLength)
  where
    (k, th) = gammaMeanVarianceToShapeScale m v

-- A variant of the log normal distribution. See Yang 2006, equation (7.23).
logNormal' :: Mean -> Variance -> Double -> Log Double
logNormal' mu var r = Exp $ negate t - e
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
uncorrelatedLogNormal :: HandleLayer -> Mean -> Variance -> PriorFunction (Tree Length a)
uncorrelatedLogNormal hl mu var = branchesWith hl (logNormal' mu var . fromLength)

-- | White noise model.
--
-- The rates are distributed according to a white noise process with given
-- variance. This is equivalent to the branch lengths measured in expected
-- number of substitutions per unit time being distributed according to a gamma
-- distribution with mean 1 and the given variance.
--
-- NOTE: The time tree has to be given because the branch length measured in
-- expected number of substitutions per unit time has to be calculated. Long
-- branches measured in time are expected to have a distribution of rates with
-- lower variance than short branches.
--
-- For example,
-- @
-- prior = whiteNoise variance timeTree rateTree
-- @
--
-- NOTE: The name white noise implies that the process is uncorrelated, but the
-- prefix is still kept to maintain consistency with the other names.
--
-- Lepage, T., Bryant, D., Philippe, H., & Lartillot, N., A general comparison
-- of relaxed molecular clock models, Molecular Biology and Evolution, 24(12),
-- 2669–2680 (2007). http://dx.doi.org/10.1093/molbev/msm193
--
-- Call 'error' if the topologies of the time and rate trees do not match.
whiteNoise :: HandleLayer -> Variance -> Tree Length a -> PriorFunction (Tree Length a)
whiteNoise hl v tTr rTr = branchesWith hl f zTr
  where
    zTr =
      fromMaybe
        (error "whiteNoise: Topologies of time and rate trees do not match.")
        (zipTrees tTr rTr)
    -- This is correct. The mean of b=tr is t, the variance of b is
    -- Var(tr) = t^2Var(r) = t^2 v/t = vt, as required in Lepage, 2006.
    f (t, r) = let k = fromLength t / v in gamma k (recip k) (fromLength r)
-- whiteNoise :: HandleLayer -> Variance -> Tree Length a -> PriorFunction (Tree Length a)
-- whiteNoise WithStem v t r =
--   -- This is correct. The mean of b=tr is t, the variance of b is
--   -- Var(tr) = t^2Var(r) = t^2 v/t = vt, as required in Lepage, 2006.
--   gamma k (recip k) (fromLength $ branch r) * whiteNoise WithoutStem v t r
--   where
--     k = fromLength (branch t) / v
-- whiteNoise WithoutStem v t r =
--   product $ zipWith (whiteNoise WithStem v) (forest t) (forest r)

-- | Auto-correlated gamma model.
--
-- Each branch length is distributed according to a gamma distribution. If the
-- parent branch exists, set the mean of the gamma distribution to the parent
-- branch length. Otherwise, use the given initial mean. The variance of the
-- gamma distribution is set to the given variance multiplied with the branch
-- length measured in unit time (see note below).
--
-- NOTE: For convenience, the mean and variance are used as parameters for this
-- relaxed molecular clock model. They are used to calculate the shape and the
-- scale of the underlying gamma distribution.
-- For example,
--
-- NOTE: The time tree has to be given because long branches are expected to
-- have a distribution of rates with higher variance than short branches. This
-- is the opposite property compared to the white noise process ('whiteNoise').
--
-- @
-- prior = autocorrelatedGamma initialMean variance timeTree rateTree
-- @
--
-- Call 'error' if the topologies of the time and rate trees do not match.
autocorrelatedGamma ::
  HandleLayer ->
  Mean ->
  Variance ->
  Tree Length a ->
  PriorFunction (Tree Length a)
autocorrelatedGamma hl mu var tTr rTr = branchesWith hl f zTr
  where
    zTr =
      fromMaybe
        (error "autocorrelatedGamma: Topologies of time and rate trees do not match.")
        (zipTrees tTr rTr)
    f (t, r) =
      let var' = fromLength t * var
          (shape, scale) = gammaMeanVarianceToShapeScale mu var'
       in gamma shape scale $ fromLength r
-- autocorrelatedGamma WithStem mu var tTr rTr =
--   gamma shape scale r * autocorrelatedGamma WithoutStem r var tTr rTr
--   where
--     t = fromLength (branch tTr)
--     r = fromLength (branch rTr)
--     var' = t * var
--     (shape, scale) = gammaMeanVarianceToShapeScale mu var'
-- autocorrelatedGamma WithoutStem mu var tTr rTr =
--   product $ zipWith (autocorrelatedGamma WithStem mu var) (forest tTr) (forest rTr)

-- | Auto-correlated log normal model.
--
-- Let \(R\) be the rate of the parent branch, and \(\mu\) and \(\sigma^2\) be
-- two given values roughly denoting mean and variance. Further, let \(t\) be
-- the length of the current branch measured in unit time. Then, the rate of the
-- current branch \(r\) is distributed according to a normal distribution with
-- mean \(\log{R} - \sigma^2 t/2) and variance \(\sigma^2t\).
--
-- The correction term is needed to ensure that the mean stays constant. See
-- Computational Molecular Evolution (Yang, 2006), Section 7.4.3.
--
-- NOTE: The time tree has to be given because long branches are expected to
-- have a distribution of rates with higher variance than short branches. This
-- is the opposite property compared to the white noise process ('whiteNoise').
--
-- For example,
-- @
-- prior = autocorrelatedLogNormal initialRate variance timeTree rateTree
-- @
--
-- Call 'error' if the topologies of the time and rate trees do not match.
autocorrelatedLogNormal ::
  HandleLayer ->
  Mean ->
  Variance ->
  Tree Length a ->
  PriorFunction (Tree Length a)
autocorrelatedLogNormal hl mu var tTr rTr = branchesWith hl f zTr
  where
    zTr =
      fromMaybe
        (error "autocorrelatedLogNormal: Topologies of time and rate trees do not match.")
        (zipTrees tTr rTr)
    f (t, r) = let var' = fromLength t * var in logNormal' mu var' $ fromLength r
-- autocorrelatedLogNormal WithStem mu var tTr rTr =
--   logNormal' mu var' r * autocorrelatedLogNormal WithoutStem r var tTr rTr
--   where
--     t = fromLength (branch tTr)
--     r = fromLength (branch rTr)
--     var' = t * var
-- autocorrelatedLogNormal WithoutStem mu var tTr rTr =
--   product $ zipWith (autocorrelatedLogNormal WithStem mu var) (forest tTr) (forest rTr)
