{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal.Hamiltonian.Common
-- Description :  Code shared by proposals based on Hamiltonian dynamics
-- Copyright   :  (c) 2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri Jun  3 10:27:12 2022.
module Mcmc.Proposal.Hamiltonian.Common
  ( -- * Basic types
    Positions,
    Momenta,
    Masses,
    LeapfrogTrajectoryLength,
    LeapfrogScalingFactor,

    -- * Tuning
    HTuneLeapfrog (..),
    HTuneMasses (..),
    HTuningConf (..),
    massesToTuningParameters,
    tuningParametersToMasses,
    tuneDiagonalMassesOnly,
    tuneAllMasses,

    -- * Structure of state
    HSpec (..),
    checkHSpecWith,

    -- * Target function
    HTarget (..),

    -- * Hamiltonian dynamics
    HMu,
    HMassesInv,
    HData (..),
    getHData,
    generateMomenta,
    exponentialKineticEnergy,

    -- * Leapfrog integrator
    Target,
    leapfrog,
  )
where

import Data.Typeable
import qualified Data.Vector as VB
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Mcmc.Jacobian
import Mcmc.Likelihood
import Mcmc.Prior
import Mcmc.Proposal
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import qualified Statistics.Covariance as S
import qualified Statistics.Function as S
import qualified Statistics.Sample as S
import System.Random.MWC

-- NOTE: Implementing the Riemannian adaptation (state-dependent mass matrix).
-- seems a little bit of an overkill.

-- | The Hamiltonian proposal acts on a vector of floating point values referred
-- to as positions.
--
-- The positions can represent the complete state or a subset of the state of
-- the Markov chain.
--
-- The length of the position vector determines the size of the squared mass
-- matrix 'Masses'.
type Positions = VS.Vector Double

-- | Internal. Momenta of the 'Positions'.
type Momenta = VS.Vector Double

-- | Parameter mass matrix.
--
-- The masses roughly describe how reluctant the particles move through the
-- state space. If a parameter has higher mass, the velocity in this direction
-- will be changed less by the provided gradient, than when the same parameter
-- has lower mass. Off-diagonal entries describe the covariance structure. If
-- two parameters are negatively correlated, their generated initial momenta are
-- likely to have opposite signs.
--
-- The matrix is square with the number of rows/columns being equal to the
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
-- \([\text{floor}(0.8L),\text{ceiling}(1.2L)]\).
--
-- For a discussion of ergodicity and reasons why randomization is important,
-- see [1] p. 15; also mentioned in [2] p. 304.
--
-- Usually set to 10, but larger values may be desirable.
--
-- NOTE: To avoid errors, the left bound of the interval has an additional hard
-- minimum of 1, and the right bound is required to be larger equal than the
-- left bound.
--
-- NOTE: Call 'error' if value is less than 1.
type LeapfrogTrajectoryLength = Int

-- | Mean of leapfrog scaling factor \(\epsilon\).
--
-- Determines the size of each leapfrog step.
--
-- To avoid problems with ergodicity, the actual leapfrog scaling factor is
-- sampled per proposal from a continuous uniform distribution over the interval
-- \((0.8\epsilon,1.2\epsilon]\).
--
-- For a discussion of ergodicity and reasons why randomization is important,
-- see [1] p. 15; also mentioned in [2] p. 304.
--
-- Usually set such that \( L \epsilon = 1.0 \), but smaller values may be
-- required if acceptance rates are low.
--
-- NOTE: Call 'error' if value is zero or negative.
type LeapfrogScalingFactor = Double

-- | Tune leapfrog parameters?
data HTuneLeapfrog
  = HNoTuneLeapfrog
  | -- | We expect that the larger the leapfrog scaling factor the lower the
    -- acceptance ratio. Consequently, if the acceptance rate is too low, the
    -- leapfrog scaling factor is decreased and vice versa. Further, the leapfrog
    -- trajectory length is scaled such that the product of the leapfrog scaling
    -- factor and leapfrog trajectory length stays roughly constant.
    HTuneLeapfrog
  deriving (Eq, Show)

-- | Tune masses?
--
-- The masses are tuned according to the (co)variances of the parameters
-- obtained from the posterior distribution over the last auto tuning interval.
data HTuneMasses
  = HNoTuneMasses
  | -- | Diagonal only: The variances of the parameters are calculated and the
    -- masses are amended using the old masses and the inverted variances. If, for
    -- a specific coordinate, the sample size is 60 or lower, or if the calculated
    -- variance is out of predefined bounds [1e-6, 1e6], the mass of the affected
    -- position is not changed.
    HTuneDiagonalMassesOnly
  | -- | All masses: The covariance matrix of the parameters is estimated and the
    -- inverted matrix (sometimes called precision matrix) is used as mass matrix.
    -- This procedure is error prone, but models with high correlations between
    -- parameters it is necessary to tune off-diagonal entries. The full mass
    -- matrix is only tuned if more than 200 samples are available. For these
    -- reasons, when tuning all masses it is recommended to use tuning settings
    -- such as
    --
    -- @
    -- BurnInWithCustomAutoTuning ([10, 20 .. 200] ++ replicate 5 500)
    -- @
    HTuneAllMasses
  deriving (Eq, Show)

-- | Tuning configuration of the Hamilton Monte Carlo proposal.
data HTuningConf = HTuningConf HTuneLeapfrog HTuneMasses
  deriving (Eq, Show)

-- If changed, also change help text of 'HTuneMasses'.
massMin :: Double
massMin = 1e-6

-- If changed, also change help text of 'HTuneMasses'.
massMax :: Double
massMax = 1e6

-- Minimal number of unique samples required for tuning the diagonal entries of
-- the mass matrix.
--
-- If changed, also change help text of 'HTuneMasses'.
samplesMinDiagonal :: Int
samplesMinDiagonal = 61

-- Minimal number of samples required for tuning all entries of the mass matrix.
--
-- If changed, also change help text of 'HTuneMasses'.
samplesMinAll :: Int
samplesMinAll = 201

getSampleSize :: VS.Vector Double -> Int
getSampleSize = VS.length . VS.uniq . S.gsort

-- Diagonal elements are variances which are strictly positive.
getNewMassDiagonalWithRescue :: Int -> Double -> Double -> Double
getNewMassDiagonalWithRescue sampleSize massOld massEstimate
  | sampleSize < samplesMinDiagonal = massOld
  -- NaN and negative masses could be errors.
  | isNaN massEstimate = massOld
  | massEstimate <= 0 = massOld
  | massMin > massNew = massMin
  | massNew > massMax = massMax
  | otherwise = massNew
  where
    massNewSqrt = recip 3 * (sqrt massOld + 2 * sqrt massEstimate)
    massNew = massNewSqrt ** 2

-- TODO @Dominik (high, issue): I do not think 'massesToTuningParameters' should
-- be exported. Something is wrong here; fix!

massesToTuningParameters :: Masses -> AuxiliaryTuningParameters
massesToTuningParameters = VB.convert . L.flatten . L.unSym

-- TODO @Dominik (high, issue): I do not think 'tuningParametersToMasses' should
-- be exported. Something is wrong here; fix!

tuningParametersToMasses ::
  -- Dimension of the mass matrix.
  Int ->
  AuxiliaryTuningParameters ->
  Masses
tuningParametersToMasses d = L.trustSym . L.reshape d . VB.convert

-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneDiagonalMassesOnly ::
  -- Number of parameters.
  Int ->
  -- Conversion from value to vector.
  (a -> Positions) ->
  -- Value vector.
  VB.Vector a ->
  -- Old masses.
  VU.Vector TuningParameter ->
  -- New masses.
  VU.Vector TuningParameter
tuneDiagonalMassesOnly dim toVec xs ts
  -- If not enough data is available, do not tune.
  | VB.length xs < samplesMinDiagonal = ts
  | otherwise =
      -- Replace the diagonal.
      massesToTuningParameters $
        L.trustSym $ massesOld - L.diag massesDiagonalOld + L.diag massesDiagonalNew
  where
    -- xs: Each vector entry contains all parameter values of one iteration.
    -- xs': Each row contains all parameter values of one iteration.
    xs' = L.fromRows $ VB.toList $ VB.map toVec xs
    sampleSizes = VS.fromList $ map getSampleSize $ L.toColumns xs'
    massesOld = L.unSym $ tuningParametersToMasses dim ts
    massesDiagonalOld = L.takeDiag massesOld
    massesDiagonalEstimate = VS.fromList $ map (recip . S.variance) $ L.toColumns xs'
    massesDiagonalNew =
      VS.zipWith3
        getNewMassDiagonalWithRescue
        sampleSizes
        massesDiagonalOld
        massesDiagonalEstimate

-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneAllMasses ::
  -- Number of parameters.
  Int ->
  -- Conversion from value to vector.
  (a -> Positions) ->
  -- Value vector.
  VB.Vector a ->
  -- Old masses.
  AuxiliaryTuningParameters ->
  -- New masses.
  AuxiliaryTuningParameters
tuneAllMasses dim toVec xs ts
  -- If not enough data is available, do not tune.
  | VB.length xs < samplesMinDiagonal = ts
  -- If not enough data is available, only the diagonal masses are tuned.
  | VB.length xs < samplesMinAll = fallbackDiagonal
  | L.rank xs' /= dim = fallbackDiagonal
  | otherwise = massesToTuningParameters $ L.trustSym massesNew
  where
    fallbackDiagonal = tuneDiagonalMassesOnly dim toVec xs ts
    -- xs: Each vector entry contains all parameter values of one iteration.
    -- xs': Each row contains all parameter values of one iteration.
    xs' = L.fromRows $ VB.toList $ VB.map toVec xs
    (_, ss, xsNormalized) = S.scale xs'
    -- sigmaNormalized = L.unSym $ either error id $ S.oracleApproximatingShrinkage xsNormalized
    sigmaNormalized = L.unSym $ either error fst $ S.graphicalLasso 0.5 xsNormalized
    sigma = S.rescaleSWith ss sigmaNormalized
    massesNew = L.inv sigma

-- | The Hamilton Monte Carlo proposal requires information about the structure
-- of the state, which is denoted as @s@.
--
-- Please also refer to the top level module documentation.
data HSpec s = HSpec
  { -- | The sample state is used for error checks and to calculate the dimension
    -- of the proposal.
    hSample :: s Double,
    -- | Extract a a subset of values to be manipulated by the Hamiltonian
    -- proposal from the complete state.
    hToVector :: s Double -> Positions,
    -- | Put those values back into the complete state.
    hFromVectorWith :: s Double -> Positions -> s Double
  }

checkHSpecWith :: Eq (s Double) => Masses -> HSpec s -> Maybe String
checkHSpecWith ms (HSpec x toVec fromVec)
  | fromVec x xVec /= x = eWith "'fromVectorWith x (toVector x) /= x' for sample state."
  | L.size xVec /= nrows = eWith "Mass matrix and 'toVector x' have different sizes for sample state."
  | otherwise = Nothing
  where
    eWith m = Just $ "checkHSpecWith: " <> m
    nrows = L.rows $ L.unSym ms
    xVec = toVec x

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

-- | Internal. Mean vector containing zeroes. We save this vector because it is
-- required when sampling from the multivariate normal distribution.
type HMu = L.Vector Double

-- | Internal. Symmetric, inverted mass matrix.
type HMassesInv = L.Herm Double

-- Internal data type containing memoized values.
data HData = HData
  { _hMu :: HMu,
    _hMassesInv :: HMassesInv
  }

-- | Compute inverted mass matrix.
--
-- Call 'error' if the determinant of the covariance matrix is negative.
getHData :: Masses -> HData
getHData ms =
  -- The multivariate normal distribution requires a positive definite matrix
  -- with positive determinant.
  if sign == 1.0
    then HData mu massesInvH
    else
      let msg =
            "getHData: Determinant of covariance matrix is negative."
       in error msg
  where
    nrows = L.rows $ L.unSym ms
    mu = L.fromList $ replicate nrows 0.0
    (massesInv, (_, sign)) = L.invlndet $ L.unSym ms
    -- NOTE: In theory we can trust that the matrix is symmetric here, because
    -- the inverse of a symmetric matrix is symmetric. However, one may want to
    -- implement a check anyways.
    massesInvH = L.trustSym massesInv

-- | Generate momenta for a new iteration.
generateMomenta ::
  HMu ->
  Masses ->
  GenIO ->
  IO Momenta
generateMomenta mu masses gen = do
  seed <- uniformM gen :: IO Int
  let momenta = L.gaussianSample seed 1 mu masses
  return $ L.flatten momenta

-- TODO (medium): Use a sparse matrix approach for the log density of the
-- multivariate normal, similar to McmcDate.

-- Exponential of kinetic energy.
exponentialKineticEnergy ::
  -- Inverted mass matrix.
  L.Herm Double ->
  -- Momenta.
  Momenta ->
  Log Double
exponentialKineticEnergy msInvH xs =
  Exp $ (-0.5) * ((xs L.<# msInv) L.<.> xs)
  where
    msInv = L.unSym msInvH

-- | Internal; Function calculating target value and gradient.
--
-- The function acts on the subset of the state manipulated by the proposal but
-- the value and gradient have to be calculated for the complete state. The
-- reason is that parameters untouched by the Hamiltonian proposal may affect
-- the result or the gradient.
--
-- Make sure that the value is calculated lazily because many times, only the
-- gradient is required.
type Target = Positions -> (Log Double, Positions)

-- | Internal; Leapfrog integrator (also used by NUTS proposal).
leapfrog ::
  Target ->
  HMassesInv ->
  --
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  --
  Positions ->
  Momenta ->
  -- | (New positions, new momenta, old target, new target).
  --
  -- Fail if state is not valid.
  Maybe (Positions, Momenta, Log Double, Log Double)
leapfrog tF hMassesInv l eps theta phi = do
  let -- The first half step of the momenta.
      (x, phiHalf) = leapfrogStepMomenta (0.5 * eps) tF theta phi
  -- L-1 full steps for positions and momenta. This gives the positions
  -- theta_{L-1}, and the momenta phi_{L-1/2}.
  (thetaLM1, phiLM1Half) <- go (l - 1) $ Just $ (theta, phiHalf)
  -- The last full step of the positions.
  let thetaL = leapfrogStepPositions hMassesInv eps thetaLM1 phiLM1Half
  -- The last half step of the momenta.
  let (x', phiL) = leapfrogStepMomenta (0.5 * eps) tF thetaL phiLM1Half
  return (thetaL, phiL, x, x')
  where
    go _ Nothing = Nothing
    go n (Just (t, p))
      | n <= 0 = Just (t, p)
      | otherwise =
          let t' = leapfrogStepPositions hMassesInv eps t p
              (x, p') = leapfrogStepMomenta eps tF t' p
           in if x /= 0.0
                then go (n - 1) $ Just $ (t', p')
                else Nothing

leapfrogStepMomenta ::
  LeapfrogScalingFactor ->
  Target ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New momenta; also return value target function to be collected at the end
  -- of the leapfrog integration.
  (Log Double, Momenta)
leapfrogStepMomenta eps tf theta phi = (x, phi + L.scale eps g)
  where
    (x, g) = tf theta

leapfrogStepPositions ::
  HMassesInv ->
  LeapfrogScalingFactor ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New positions.
  Positions
leapfrogStepPositions hMassesInv eps theta phi = theta + (L.unSym hMassesInvEps L.#> phi)
  where
    hMassesInvEps = L.scale eps hMassesInv
