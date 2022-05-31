{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}

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
-- The Hamiltonian Monte Carlo (HMC) proposal.
--
-- For references, see:
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
-- Notes on implementation:
--
-- TODO (high): Clean up!
--
-- - The 'Gradient' of the log posterior needs to be provided. Like so, the user
--   can use automatic or manual differentiation, depending on the problem at
--   hand.
--
--   Further, the gradient of the complete state needs to be provided, even
--   though the proposal may only act on a sub-set of the complete state. In
--   particular, do not use 'liftProposalWith', 'liftProposal', or '(@~)' with
--   the Hamiltonian Monte Carlo proposal; instead, see the documentation of
--   'HSpec'.
--
-- - The Hamiltonian proposal acts on a vector of storable 'Positions'.
--   Functions converting the state to and from this vector have to be provided.
--   See 'HSpec'.
--
-- - The desired acceptance rate is 0.65, although the dimension of the proposal
--   is high.
--
-- - The speed of this proposal changes drastically with the leapfrog trajectory
--   length and the leapfrog scaling factor. Hence, the speed will change during
--   burn in.
module Mcmc.Proposal.Hamiltonian
  ( Positions,
    Momenta,
    Masses,
    LeapfrogTrajectoryLength,
    LeapfrogScalingFactor,
    HTuneLeapfrog (..),
    HTuneMasses (..),
    HTuningConf (..),
    HTuningSpec,
    hTuningSpec,
    HSpec (..),
    HTarget (..),
    HMassesInv,
    Target,
    leapfrog,
    hamiltonian,
  )
where

import Data.Bifunctor
import Data.Typeable
import qualified Data.Vector as VB
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Mcmc.Jacobian
import Mcmc.Likelihood
import Mcmc.Prior
import Mcmc.Proposal
import Numeric.AD.Double
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import Numeric.MathFunctions.Constants
import qualified Statistics.Covariance as S
import qualified Statistics.Function as S
import qualified Statistics.Sample as S
import System.Random.MWC

-- TODO (high): Clean up!

-- NOTE: Implementing the Riemannian adaptation (state-dependent mass matrix).
-- seems a little bit of an overkill.

-- | The Hamiltonian proposal acts on a vector of floating point values referred
-- to as positions.
--
-- The positions can represent the complete state or a subset of the state of
-- the Markov chain.
type Positions = L.Vector Double

-- | Internal. Momenta of the 'Positions'.
type Momenta = L.Vector Double

-- | Parameter mass matrix.
--
-- The masses roughly describe how reluctant the particles move through the
-- state space. If a parameter has higher mass, the velocity in this direction
-- will be changed less by the provided gradient, than when the same parameter
-- has lower mass. Off-diagonal entries describe the covariance structure. If
-- two parameters are negatively correlated, their generated initial momenta are
-- likely to have opposite signs.
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

-- | Complete tuning specification of the Hamilton Monte Carlo proposal.
--
-- Includes tuning parameters and tuning configuration.
data HTuningSpec = HTuningSpec
  { hMasses :: Masses,
    hLeapfrogTrajectoryLength :: LeapfrogTrajectoryLength,
    hLeapfrogScalingFactor :: LeapfrogScalingFactor,
    hTuningConf :: HTuningConf
  }
  deriving (Show)

-- | See 'HTuningSpec'.
--
-- Return 'Left' if an error is found.
hTuningSpec ::
  Masses ->
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  HTuningConf ->
  Either String HTuningSpec
hTuningSpec masses l eps c
  | any (<= 0) diagonalMasses = eWith "Some diagonal entries of the mass matrix are zero or negative."
  | nrows /= ncols = eWith "Mass matrix is not square."
  | l < 1 = eWith "Leapfrog trajectory length is zero or negative."
  | eps <= 0 = eWith "Leapfrog scaling factor is zero or negative."
  | otherwise = Right $ HTuningSpec masses l eps c
  where
    eWith m = Left $ "hTuningSpec: " <> m
    ms = L.unSym masses
    diagonalMasses = L.toList $ L.takeDiag ms
    nrows = L.rows ms
    ncols = L.cols ms

-- | The Hamilton Monte Carlo proposal requires information about the structure
-- of the state, which is denoted as @s@.
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

-- | The target is composed of the prior, likelihood, and jacobian functions.
--
-- The structure of the state is denoted as @s@.
data HTarget s = HTarget
  { -- | Function computing the log prior.
    hPrior :: forall a. (RealFloat a, Typeable a) => PriorFunctionG (s a) a,
    -- | Function computing the log likelihood.
    hLikelihood :: forall a. (RealFloat a, Typeable a) => LikelihoodFunctionG (s a) a,
    -- | Function computing the log of the Jacobian.
    hJacobian :: forall a. (RealFloat a, Typeable a) => Maybe (JacobianFunctionG (s a) a)
  }

checkHSpecWith :: Eq (s Double) => HTuningSpec -> HSpec s -> Maybe String
checkHSpecWith tspec (HSpec x toVec fromVec)
  | fromVec x xVec /= x = eWith "'fromVectorWith x (toVector x) /= x' for sample state."
  | L.size xVec /= nrows = eWith "Mass matrix and 'toVector x' have different sizes for sample state."
  | otherwise = Nothing
  where
    eWith m = Just $ "checkHSpecWith: " <> m
    ms = L.unSym $ hMasses tspec
    nrows = L.rows ms
    xVec = toVec x

-- Internal. Mean vector containing zeroes.
type HMu = L.Vector Double

-- | Internal. Symmetric, inverted mass matrix.
type HMassesInv = L.Herm Double

-- Internal. Logarithm of the determinant of the mass matrix.
type HLogDetMasses = Double

-- Internal data type containing memoized values.
data HData = HData
  { _hMu :: HMu,
    _hMassesInv :: HMassesInv,
    _hLogDetMasses :: HLogDetMasses
  }

-- Call 'error' if the determinant of the covariance matrix is negative.
getHData :: HTuningSpec -> HData
getHData s =
  -- The multivariate normal distribution requires a positive definite matrix
  -- with positive determinant.
  if sign == 1.0
    then HData mu massesInvH logDetMasses
    else
      let msg =
            "hamiltonianSimple: Determinant of covariance matrix is negative."
              <> " The logarithm of the absolute value of the determinant is: "
              <> show logDetMasses
              <> "."
       in error msg
  where
    ms = hMasses s
    nrows = L.rows $ L.unSym ms
    mu = L.fromList $ replicate nrows 0.0
    (massesInv, (logDetMasses, sign)) = L.invlndet $ L.unSym ms
    -- In theory we can trust that the matrix is symmetric here, because the
    -- inverse of a symmetric matrix is symmetric. However, one may want to
    -- implement a check anyways.
    massesInvH = L.trustSym massesInv

generateMomenta ::
  -- Provided so that it does not have to be recreated.
  HMu ->
  Masses ->
  GenIO ->
  IO Momenta
generateMomenta mu masses gen = do
  seed <- uniformM gen :: IO Int
  let momenta = L.gaussianSample seed 1 mu masses
  return $ L.flatten momenta

-- Prior distribution of momenta.
--
-- Log of density of multivariate normal distribution with given parameters.
-- https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Density_function.
logDensityMultivariateNormal ::
  -- Mean vector.
  L.Vector Double ->
  -- Inverted covariance matrix.
  L.Herm Double ->
  -- Logarithm of the determinant of the covariance matrix.
  Double ->
  -- Value vector.
  L.Vector Double ->
  Log Double
logDensityMultivariateNormal mu sigmaInvH logDetSigma xs =
  Exp $ c + (-0.5) * (logDetSigma + ((dxs L.<# sigmaInv) L.<.> dxs))
  where
    dxs = xs - mu
    k = fromIntegral $ L.size mu
    c = negate $ m_ln_sqrt_2_pi * k
    sigmaInv = L.unSym sigmaInvH

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
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  Positions ->
  Momenta ->
  -- | Fail if state is not valid.
  Maybe (Positions, Momenta, Log Double)
leapfrog tF hMassesInv l eps theta phi = do
  let -- The first half step of the momenta.
      (_, phiHalf) = leapfrogStepMomenta (0.5 * eps) tF theta phi
  -- L-1 full steps for positions and momenta. This gives the positions
  -- theta_{L-1}, and the momenta phi_{L-1/2}.
  (thetaLM1, phiLM1Half) <- go (l - 1) $ Just $ (theta, phiHalf)
  -- The last full step of the positions.
  let thetaL = leapfrogStepPositions hMassesInv eps thetaLM1 phiLM1Half
  -- The last half step of the momenta.
  let (x, phiL) = leapfrogStepMomenta (0.5 * eps) tF thetaL phiLM1Half
  return (thetaL, phiL, x)
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

massesToTuningParameters :: Masses -> AuxiliaryTuningParameters
massesToTuningParameters = VB.convert . L.flatten . L.unSym

tuningParametersToMasses ::
  -- Dimension of the mass matrix.
  Int ->
  AuxiliaryTuningParameters ->
  Masses
tuningParametersToMasses d = L.trustSym . L.reshape d . VB.convert

hTuningParametersToTuningSpec ::
  HTuningSpec ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String HTuningSpec
hTuningParametersToTuningSpec (HTuningSpec ms l e c) t ts
  | not nTsOK = Left "hTuningParametersToSettings: Auxiliary variables dimension mismatch."
  | otherwise = Right $ HTuningSpec msTuned lTuned eTuned c
  where
    d = L.rows $ L.unSym ms
    (HTuningConf tlf tms) = c
    nTsOK =
      let nTs = VU.length ts
       in case tms of
            HNoTuneMasses -> nTs == 0
            _ -> nTs == d * d
    msTuned = case tms of
      HNoTuneMasses -> ms
      _ -> tuningParametersToMasses d ts
    -- The larger epsilon, the larger the proposal step size and the lower the
    -- expected acceptance ratio.
    --
    -- Further, we roughly keep \( L * \epsilon = 1.0 \). The equation is not
    -- correct, because we pull L closer to the original value to keep the
    -- runtime somewhat acceptable.
    (lTuned, eTuned) = case tlf of
      HNoTuneLeapfrog -> (l, e)
      HTuneLeapfrog -> (ceiling $ fromIntegral l / (t ** 0.9) :: Int, t * e)

hamiltonianSimpleWithTuningParameters ::
  Traversable s =>
  HTuningSpec ->
  HSpec s ->
  (s Double -> Target) ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (ProposalSimple (s Double))
hamiltonianSimpleWithTuningParameters tspec hspec targetWith t ts = do
  tspec' <- hTuningParametersToTuningSpec tspec t ts
  pure $ hamiltonianSimple tspec' hspec targetWith

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianSimple'.
hamiltonianSimpleWithMemoizedCovariance ::
  Traversable s =>
  HTuningSpec ->
  HSpec s ->
  HData ->
  (s Double -> Target) ->
  ProposalSimple (s Double)
hamiltonianSimpleWithMemoizedCovariance tspec hspec dt targetWith x g = do
  phi <- generateMomenta mu masses g
  lRan <- uniformR (lL, lR) g
  eRan <- uniformR (eL, eR) g
  case leapfrog (targetWith x) massesInv lRan eRan theta phi of
    Nothing -> return (x, 0.0, 1.0)
    Just (theta', phi', _) ->
      let -- Prior of momenta.
          prPhi = logDensityMultivariateNormal mu massesInv logDetMasses phi
          prPhi' = logDensityMultivariateNormal mu massesInv logDetMasses phi'
          kernelR = prPhi' / prPhi
       in -- NOTE: For example, Neal page 12: In order for the Hamiltonian proposal
          -- to be in detailed balance, the momenta have to be negated before
          -- proposing the new value. That is, the negated momenta would guide the
          -- chain back to the previous state. However, we are only interested in
          -- the positions, and are not even storing the momenta.
          return (fromVec x theta', kernelR, 1.0)
  where
    (HTuningSpec masses l e _) = tspec
    (HSpec _ toVec fromVec) = hspec
    theta = toVec x
    lL = maximum [1 :: Int, floor $ (0.8 :: Double) * fromIntegral l]
    lR = maximum [lL, ceiling $ (1.2 :: Double) * fromIntegral l]
    eL = 0.8 * e
    eR = 1.2 * e
    (HData mu massesInv logDetMasses) = dt

hamiltonianSimple ::
  Traversable s =>
  HTuningSpec ->
  HSpec s ->
  (s Double -> Target) ->
  ProposalSimple (s Double)
hamiltonianSimple tspec hspec = hamiltonianSimpleWithMemoizedCovariance tspec hspec (getHData tspec)

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

-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneDiagonalMassesOnly ::
  Int ->
  (a -> Positions) ->
  AuxiliaryTuningFunction a
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
  Int ->
  (a -> Positions) ->
  AuxiliaryTuningFunction a
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

-- | Hamiltonian Monte Carlo proposal.
hamiltonian ::
  (Eq (s Double), Traversable s) =>
  HTuningSpec ->
  HSpec s ->
  HTarget s ->
  PName ->
  PWeight ->
  Proposal (s Double)
hamiltonian tspec hspec htarget n w = case checkHSpecWith tspec hspec of
  Just err -> error err
  Nothing ->
    let -- Misc.
        desc = PDescription "Hamiltonian Monte Carlo (HMC)"
        (HSpec sample toVec fromVec) = hspec
        (HTarget prF lhF mJF) = htarget
        dim = (L.size $ toVec sample)
        pDim = PSpecial dim 0.65
        -- Vectorize and derive the target function.
        tF y = case mJF of
          Nothing -> prF y * lhF y
          Just jF -> prF y * lhF y * jF y
        tFnG = grad' (ln . tF)
        targetWith x = bimap Exp toVec . tFnG . fromVec x
        ps = hamiltonianSimple tspec hspec targetWith
        hamiltonianWith = Proposal n desc PSlow pDim w ps
        -- Tuning.
        ts = massesToTuningParameters $ hMasses tspec
        tSet@(HTuningConf tlf tms) = hTuningConf tspec
        tFun = case tlf of
          HNoTuneLeapfrog -> noTuningFunction
          HTuneLeapfrog -> defaultTuningFunctionWith pDim
        tFunAux = case tms of
          HNoTuneMasses -> noAuxiliaryTuningFunction
          HTuneDiagonalMassesOnly -> tuneDiagonalMassesOnly dim toVec
          HTuneAllMasses -> tuneAllMasses dim toVec
     in case tSet of
          (HTuningConf HNoTuneLeapfrog HNoTuneMasses) -> hamiltonianWith Nothing
          _ ->
            let tuner = Tuner 1.0 tFun ts tFunAux (hamiltonianSimpleWithTuningParameters tspec hspec targetWith)
             in hamiltonianWith $ Just tuner
