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
-- The Hamiltonian Monte Carlo ('HMC') proposal.
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
  ( HMC (..),
    hmc,
  )
where

import Mcmc.Prior
import Mcmc.Proposal
import System.Random.MWC

-- | Specifications for Hamilton Monte Carlo proposal.
data HMC f = HMC
  { -- | Gradient of the posterior function.
    hmcGradient :: f Double -> f Double,
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
    hmcMasses :: f Double,
    -- | Number \(L\) of leapfrog steps per proposal.
    --
    -- Usually set to 10.
    hmcNLeapFrogSteps :: Int,
    -- | Scaling factor \(\epsilon\) of the leapfrog steps.
    --
    -- Usually set such that \( L \epsilon = 1.0 \).
    hmcLeapFrogScalingFactor :: Double
  }

generateMomenta ::
  -- Masses.
  f Double ->
  GenIO ->
  IO (f Double)
generateMomenta = undefined

priorMomenta ::
  f Double ->
  f Double ->
  Prior
priorMomenta = undefined

leapfrog ::
  Applicative f =>
  -- Number of leapfrog steps \(L\).
  Int ->
  -- Scaling factor \(\epsilon\).
  Double ->
  -- Masses.
  f Double ->
  -- Gradient.
  (f Double -> f Double) ->
  -- Positions.
  f Double ->
  -- Momenta.
  f Double ->
  -- (Positions', Momenta').
  (f Double, f Double)
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
  -- Scaling factor \(\epsilon\).
  Double ->
  -- Gradient.
  (f Double -> f Double) ->
  -- Current positions.
  f Double ->
  -- Current momenta.
  f Double ->
  -- New momenta.
  f Double
leapfrogStepMomenta xi eps grad theta phi = phi .+. ((xi * eps) .* grad theta)

leapfrogStepPositions ::
  Applicative f =>
  -- Scaling factor \(\epsilon\).
  Double ->
  -- Masses.
  f Double ->
  -- Current positions.
  f Double ->
  -- Current momenta.
  f Double ->
  f Double
leapfrogStepPositions eps masses theta phi = theta .+. (mReversedScaled .*. phi)
  where
    mReversedScaled = (* eps) . (** (-1)) <$> masses

-- Element-wise vector-vector addition.
(.+.) :: Applicative f => f Double -> f Double -> f Double
(.+.) xs ys = (+) <$> xs <*> ys

-- Element-wise vector-vector multiplication.
(.*.) :: Applicative f => f Double -> f Double -> f Double
(.*.) xs ys = (*) <$> xs <*> ys

-- Scalar-vector multiplication.
(.*) :: Applicative f => Double -> f Double -> f Double
(.*) x ys = (* x) <$> ys

-- phi half update
-- (l-1) theta and phi full updates
-- phi half update

hmcSimpleWith ::
  Applicative f =>
  HMC f ->
  TuningParameter ->
  f Double ->
  GenIO ->
  IO (f Double, KernelRatio, Jacobian)
hmcSimpleWith s t theta g = do
  phi <- generateMomenta masses g
  let (theta', phi') = leapfrog lTuned epsTuned masses gradient theta phi
      prPhi = priorMomenta masses phi
      prPhi' = priorMomenta masses phi'
      kernelR = prPhi' / prPhi
  return (theta', kernelR, 1.0)
  where
    l = hmcNLeapFrogSteps s
    -- The larger epsilon, the larger the proposal step size and the lower the
    -- acceptance ratio.
    --
    -- Further, lTuned * epsTuned should be approximately constant.
    lTuned = ceiling $ fromIntegral l / t
    epsTuned = t * hmcLeapFrogScalingFactor s
    masses = hmcMasses s
    gradient = hmcGradient s

-- | Hamiltonian Monte Carlo proposal.
hmc ::
  (Foldable f, Applicative f) =>
  -- | The sample state is used to calculate the dimension of the proposal.
  f Double ->
  HMC f ->
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
