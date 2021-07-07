{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      :  Mcmc.Proposal.Hamiltonian
-- Description :  Hamiltonian Monte Carlo proposal
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Mon Jul  5 12:59:42 2021.
--
-- The Hamiltonian Monte Carlo (HMC, see 'hmc') proposal.
--
-- For references, see:
--
-- - Chapter 5 of Handbook of Monte Carlo: Neal, R. M., MCMC Using Hamiltonian
--   Dynamics, In S. Brooks, A. Gelman, G. Jones, & X. Meng (Eds.), Handbook of
--   Markov Chain Monte Carlo (2011), CRC press.
--
-- - Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B., Bayesian data
--   analysis (2014), CRC Press.
--
-- - Review by Betancourt and notes: Betancourt, M., A conceptual introduction
--   to Hamiltonian Monte Carlo, arXiv, 1701–02434 (2017).
--
-- NOTE: This module is in development.
module Mcmc.Proposal.Hamiltonian
  ( Gradient,
    Masses,
    LeapfrogTrajectoryLength,
    LeapfrogScalingFactor,
    HmcSettings (..),
    hmc,
  )
where

import Data.Foldable
import Mcmc.Prior
import Mcmc.Proposal
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC

-- TODO: Allow random leapfrog trajectory lengths and leapfrog scaling factors
-- (sample one l and one epsilon per proposal, not per leapfrog step).

-- | Gradient of the posterior function.
--
-- NOTE: Not the log gradient in log domain.
type Gradient f = f Double -> f (Log Double)

-- | Masses.
--
-- The masses describe how fast the sampler moves in each direction of the
-- state space. High masses result in low momenta, and vice versa. The
-- proposal is more efficient if masses are assigned according to the
-- (co)-variance structure of the posterior. However, the proposal works
-- fine when all masses are set to 1.0.
--
-- NOTE: At the moment, it is impossible to set off-diagonal elements of the
-- mass matrix.
type Masses f = f Double

-- | Leapfrog trajectory length \(L\).
--
-- Number of leapfrog steps per proposal.
--
-- Usually set to 10.
type LeapfrogTrajectoryLength = Int

-- | Leapfrog scaling factor \(\epsilon\).
--
-- Determines the size of each leapfrog step.
--
-- Usually set such that \( L \epsilon = 1.0 \).
type LeapfrogScalingFactor = Double

-- State containing parameters.
type Positions f = f Double

-- Momenta or velocities of the parameters.
--
-- TODO: Check if these are actually the momenta, or the velocities (without
-- mass component).
type Momenta f = f Double

-- | Specifications for Hamilton Monte Carlo proposal.
data HmcSettings f = HmcSettings
  { hmcGradient :: Gradient f,
    hmcMasses :: Masses f,
    hmcLeapfrogTrajectoryLength :: LeapfrogTrajectoryLength,
    hmcLeapfrogScalingFactor :: LeapfrogScalingFactor
  }

generateMomenta ::
  Traversable f =>
  Masses f ->
  GenIO ->
  IO (Momenta f)
generateMomenta masses gen = traverse (generateWith gen) masses
  where
    generateWith g m = let d = normalDistr 0 m in genContVar d g

priorMomenta ::
  (Applicative f, Foldable f) =>
  Masses f ->
  Momenta f ->
  Prior
priorMomenta masses phi = foldl' (*) 1.0 $ f <$> masses <*> phi
  where
    f m p = let d = normalDistr 0 m in Exp $ logDensity d p

leapfrog ::
  Applicative f =>
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  Masses f ->
  Gradient f ->
  Positions f ->
  Momenta f ->
  -- (Positions', Momenta').
  (Positions f, Momenta f)
leapfrog l eps masses grad theta phi = (thetaL, phiL)
  where
    -- The first half step of the momenta.
    phiHalf = leapfrogStepMomenta 0.5 eps grad theta phi
    -- L-1 full steps. This gives the positions theta_{L-1}, and the momenta
    -- phi_{L-1/2}.
    (thetaLM1, phiLM1Half) = go (l - 1) (theta, phiHalf)
    -- The last full step of the positions.
    thetaL = leapfrogStepPositions eps masses thetaLM1 phiLM1Half
    -- The last half step of the momenta.
    phiL = leapfrogStepMomenta 0.5 eps grad thetaL phiLM1Half
    go 0 (t, p) = (t, p)
    go n (t, p) =
      let t' = leapfrogStepPositions eps masses t p
       in go (n - 1) (t', leapfrogStepMomenta 1.0 eps grad t' p)

leapfrogStepMomenta ::
  Applicative f =>
  -- Size of step (half or full step).
  Double ->
  LeapfrogScalingFactor ->
  Gradient f ->
  -- Current positions.
  Positions f ->
  -- Current momenta.
  Momenta f ->
  -- New momenta.
  Momenta f
leapfrogStepMomenta xi eps grad theta phi = phi .+. ((xi * eps) .* grad theta)
  where
    -- Scalar-vector multiplication.
    (.*) :: Applicative f => Double -> f (Log Double) -> f Double
    (.*) x ys = (* x) . ln <$> ys

leapfrogStepPositions ::
  Applicative f =>
  LeapfrogScalingFactor ->
  Masses f ->
  -- Current positions.
  Positions f ->
  -- Current momenta.
  Momenta f ->
  Positions f
leapfrogStepPositions eps masses theta phi = theta .+. (mReversedScaled .*. phi)
  where
    mReversedScaled = (* eps) . (** (-1)) <$> masses

-- Element-wise vector-vector addition.
(.+.) :: Applicative f => f Double -> f Double -> f Double
(.+.) xs ys = (+) <$> xs <*> ys

-- Element-wise vector-vector multiplication.
(.*.) :: Applicative f => f Double -> f Double -> f Double
(.*.) xs ys = (*) <$> xs <*> ys

-- phi half update
-- (l-1) theta and phi full updates
-- phi half update

hmcSimpleWith ::
  (Applicative f, Traversable f, Show (f Double)) =>
  HmcSettings f ->
  TuningParameter ->
  Positions f ->
  GenIO ->
  IO (Positions f, KernelRatio, Jacobian)
hmcSimpleWith s t theta g = do
  phi <- generateMomenta masses g
  let (theta', phi') = leapfrog lTuned epsTuned masses gradient theta phi
      prPhi = priorMomenta masses phi
      -- XXX: In Neal page 12, it is stated that the momenta have to be negated
      -- before proposing the new value to make the proposal symmetrical. This
      -- does not change anything here since the prior involves normal
      -- distributions centered around 0. However, if the multivariate normal
      -- distribution is used, it makes a difference.
      prPhi' = priorMomenta masses phi'
      kernelR = prPhi' / prPhi
  return (theta', kernelR, 1.0)
  where
    l = hmcLeapfrogTrajectoryLength s
    -- The larger epsilon, the larger the proposal step size and the lower the
    -- acceptance ratio.
    --
    -- Further, lTuned * epsTuned should be approximately constant.
    --
    -- TODO: Improve tuning. Leaving l*eps constant leads to very large l and
    -- very slow proposals.
    --
    -- lTuned = ceiling $ fromIntegral l / t
    lTuned = l
    epsTuned = t * hmcLeapfrogScalingFactor s
    masses = hmcMasses s
    gradient = hmcGradient s

-- | Hamiltonian Monte Carlo proposal.
hmc ::
  (Applicative f, Traversable f, Show (f Double)) =>
  -- | The sample state is used to calculate the dimension of the proposal.
  f Double ->
  HmcSettings f ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (f Double)
hmc x s = createProposal hmcDescription (hmcSimpleWith s) d
  where
    hmcDescription = PDescription "Hamiltonian Monte Carlo (HMC)"
    -- TODO: The dimension is correct, but the calculated optimal acceptance
    -- rate differs from the actual optimal rate, which is 65%).
    d = PDimension $ length x
