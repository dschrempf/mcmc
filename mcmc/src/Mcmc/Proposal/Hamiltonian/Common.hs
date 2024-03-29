{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal.Hamiltonian.Common
-- Description :  Code shared by proposals based on Hamiltonian dynamics
-- Copyright   :  2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri Jun  3 10:27:12 2022.
--
-- See "Mcmc.Proposal.Hamiltonian.Hamiltonian".
--
-- References:
--
-- - [1] Chapter 5 of Handbook of Monte Carlo: Neal, R. M., MCMC Using
--   Hamiltonian Dynamics, In S. Brooks, A. Gelman, G. Jones, & X. Meng (Eds.),
--   Handbook of Markov Chain Monte Carlo (2011), CRC press.
--
-- - [2] Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B., Bayesian data
--   analysis (2014), CRC Press.
--
-- - [3] Review by Betancourt and notes: Betancourt, M., A conceptual
--   introduction to Hamiltonian Monte Carlo, arXiv, 1701–02434 (2017).
--
-- - [4] Matthew D. Hoffman, Andrew Gelman (2014) The No-U-Turn Sampler:
--   Adaptively Setting Path Lengths in Hamiltonian Monte Carlo, Journal of
--   Machine Learning Research.
module Mcmc.Proposal.Hamiltonian.Common
  ( -- * Basic types
    Positions,
    Momenta,
    Masses,
    LeapfrogTrajectoryLength,
    LeapfrogScalingFactor,
    LeapfrogSimulationLength,

    -- * Tuning
    HTuneLeapfrog (..),
    HTuneMasses (..),
    HTuningConf (..),

    -- * Structure of state
    HStructure (..),

    -- * Target function
    HTarget (..),
  )
where

import Data.Typeable
import qualified Data.Vector.Storable as VS
import Mcmc.Jacobian
import Mcmc.Likelihood
import Mcmc.Prior
import qualified Numeric.LinearAlgebra as L

-- NOTE: Implementing the Riemannian adaptation (state-dependent mass matrix).
-- seems a little bit of an overkill. See also
-- https://discourse.mc-stan.org/t/riemann-manifold-hmc-in-stan/19466/5.

-- | The Hamiltonian proposal acts on a vector of floating point values referred
-- to as positions.
--
-- The positions can represent the complete state or a subset of the state of
-- the Markov chain; see 'HStructure'.
--
-- The length of the position vector determines the size of the squared mass
-- matrix; see 'Masses'.
type Positions = VS.Vector Double

-- | Momenta of the 'Positions'.
type Momenta = VS.Vector Double

-- | Mass matrix.
--
-- The masses roughly describe how reluctant the particles move through the
-- state space. If a parameter has higher mass the velocity in this direction
-- will be changed less by the provided gradient than when the same parameter
-- has lower mass. Off-diagonal entries describe the covariance structure. If
-- two parameters are negatively correlated, their generated initial momenta are
-- likely to have opposite signs.
--
-- The matrix is square with the number of rows and columns being equal to the
-- length of 'Positions'.
--
-- The proposal is more efficient if masses are assigned according to the
-- inverse (co)-variance structure of the posterior function. That is,
-- parameters changing on larger scales should have lower masses than parameters
-- changing on lower scales. In particular, the optimal entries of the diagonal
-- of the mass matrix are the inverted variances of the parameters distributed
-- according to the posterior function.
--
-- Of course, the scales of the parameters of the posterior function are usually
-- unknown. Often, it is sufficient to
--
-- - set the diagonal entries of the mass matrix to identical values roughly
--   scaled with the inverted estimated average variance of the posterior
--   function; or even to
--
-- - set all diagonal entries of the mass matrix to 1.0, and all other entries
--   to 0.0, and trust the tuning algorithm (see 'HTuningConf') to find the
--   correct values.
type Masses = L.Herm Double

-- | Mean leapfrog trajectory length \(L\).
--
-- Number of leapfrog steps per proposal.
--
-- To avoid problems with ergodicity, the actual number of leapfrog steps is
-- sampled per proposal from a discrete uniform distribution over the interval
-- \([\text{floor}(0.9L),\text{ceiling}(1.1L)]\).
--
-- For a discussion of ergodicity and reasons why randomization is important,
-- see [1] p. 15; also mentioned in [2] p. 304.
--
-- To avoid errors, the left bound of the interval has an additional hard
-- minimum of 1, and the right bound is required to be larger equal than the
-- left bound.
type LeapfrogTrajectoryLength = Int

-- | Mean of leapfrog scaling factor \(\epsilon\).
--
-- Determines the size of each leapfrog step.
--
-- To avoid problems with ergodicity, the actual leapfrog scaling factor is
-- sampled per proposal from a continuous uniform distribution over the interval
-- \((0.9\epsilon,1.1\epsilon]\).
--
-- For a discussion of ergodicity and reasons why randomization is important,
-- see [1] p. 15; also mentioned in [2] p. 304.
--
-- Call 'error' if value is zero or negative.
type LeapfrogScalingFactor = Double

-- | Product of 'LeapfrogTrajectoryLength' and 'LeapfrogScalingFactor'.
--
-- A good value is hard to find and varies between applications. For example,
-- see Figure 6 in [4].
--
-- Call 'error' if value is zero or negative.
type LeapfrogSimulationLength = Double

-- | Tune leapfrog parameters?
data HTuneLeapfrog
  = HNoTuneLeapfrog
  | -- | We expect that the larger the leapfrog scaling factor the lower the
    -- acceptance ratio. Consequently, if the acceptance rate is too low, the
    -- leapfrog scaling factor is decreased and vice versa. Further, the leapfrog
    -- trajectory length is scaled such that the product of the leapfrog scaling
    -- factor and leapfrog trajectory length stays constant.
    HTuneLeapfrog
  deriving (Eq, Show)

-- | Tune masses?
--
-- The masses are tuned according to the (co)variances of the parameters
-- obtained from the posterior distribution over the last auto tuning interval.
data HTuneMasses
  = HNoTuneMasses
  | -- | Diagonal only: The variances of the parameters are calculated and the
    -- masses are amended using the old masses and the inverted variances. If,
    -- for a specific coordinate, the sample size is 60 or lower, or if the
    -- calculated variance is out of predefined bounds [1e-8, 1e8], the mass of
    -- the affected position is not changed.
    HTuneDiagonalMassesOnly
  | -- | All masses: The covariance matrix of the parameters is estimated and
    -- the inverted matrix (sometimes called precision matrix) is used as mass
    -- matrix. This procedure is error prone, but models with high correlations
    -- between parameters strongly profit from tuning off-diagonal entries. The
    -- full mass matrix is only tuned if equal or more than "(n_masses +
    -- max(n_masses, 61)" samples are available. For these reasons, when tuning
    -- all masses it is recommended to use tuning settings such as
    --
    -- @
    -- BurnInWithCustomAutoTuning [10, 20 .. 200] [200, 220 .. 500]
    -- @
    HTuneAllMasses
  deriving (Eq, Show)

-- | Tuning configuration of the Hamilton Monte Carlo proposal.
data HTuningConf = HTuningConf HTuneLeapfrog HTuneMasses
  deriving (Eq, Show)

-- | The Hamilton Monte Carlo proposal requires information about the structure
-- of the state, which is denoted as @s@.
--
-- Please also refer to the top level module documentation of
-- "Mcmc.Proposal.Hamiltonian.Hamiltonian".
data HStructure s = HStructure
  { -- | The sample state is used for error checks and to calculate the dimension
    -- of the proposal.
    hSample :: s Double,
    -- | Extract a a subset of values to be manipulated by the Hamiltonian
    -- proposal from the complete state.
    hToVector :: s Double -> Positions,
    -- | Put those values back into the complete state.
    hFromVectorWith :: s Double -> Positions -> s Double
  }

-- | The target is composed of the prior, likelihood, and jacobian functions.
--
-- The structure of the state is denoted as @s@.
--
-- Please also refer to the top level module documentation.
data HTarget s = HTarget
  { -- | Function computing the log prior.
    hPrior :: forall a. (RealFloat a, Typeable a) => Maybe (PriorFunctionG (s a) a),
    -- | Function computing the log likelihood.
    hLikelihood :: forall a. (RealFloat a, Typeable a) => LikelihoodFunctionG (s a) a,
    -- | Function computing the log of the Jacobian.
    hJacobian :: forall a. (RealFloat a, Typeable a) => Maybe (JacobianFunctionG (s a) a)
  }
