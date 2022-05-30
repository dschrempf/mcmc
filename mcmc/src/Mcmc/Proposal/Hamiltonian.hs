{-# LANGUAGE TupleSections #-}
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
-- NOTE on implementation:
--
-- - The implementation assumes the existence of the 'Gradient'. Like so, the
--   user can use automatic or manual differentiation, depending on the problem
--   at hand.
--
-- - The Hamiltonian proposal acts on a vector of storable 'Position' variables.
--   Functions converting the state to and from this vector have to be provided.
--   See 'HSpec'.
--
-- - The desired acceptance rate is 0.65, although the dimension of the proposal
--   is high.
--
-- - The speed of this proposal can change drastically when tuned because the
--   leapfrog trajectory length is changed.
--
-- - The Hamiltonian proposal is agnostic of the actual prior and likelihood
--   functions, and so, points with zero posterior probability cannot be
--   detected. This affects models with constrained parameters. See Gelman p.
--   303. This problem can be ameliorated by providing a 'Validate' function so
--   that the proposal can gracefully fail as soon as the state becomes invalid.
module Mcmc.Proposal.Hamiltonian
  ( Positions,
    Momenta,
    Gradient,
    Validate,
    Masses,
    LeapfrogTrajectoryLength,
    LeapfrogScalingFactor,
    HTuneLeapfrog (..),
    HTuneMasses (..),
    HTuningConf (..),
    HTuningSpec,
    hTuningSpec,
    HSpec (..),
    leapfrog,
    hamiltonian,
  )
where

import qualified Data.Vector as VB
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Mcmc.Proposal
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import Numeric.MathFunctions.Constants
import qualified Statistics.Covariance as S
import qualified Statistics.Function as S
import qualified Statistics.Sample as S
import System.Random.MWC

-- TODO (high): The Hamiltonian proposals needs to know about the Jacobian. This
-- is because the Jacobian has to be part of the gradient (at least I think so;
-- check this!). Only then, the Hamiltonian proposal can be combined with other
-- proposals in a setting with a non-unit Jacobian.

-- NOTE: Implementing the Riemannian adaptation (state-dependent mass matrix).
-- seems a little bit of an overkill.

-- | The Hamiltonian proposal acts on a vector of floating point values referred
-- to as positions.
--
-- The positions can represent the complete state or a subset of the state of
-- the Markov chain.
type Positions = L.Vector Double

-- Internal. Momenta of the 'Positions'.
type Momenta = L.Vector Double

-- | Gradient of the log posterior function.
--
-- The gradient has to be provided for the complete state. The reason is that
-- the gradient may change if parameters untouched by the Hamiltonian proposal
-- are altered by other proposals.
type Gradient a = a -> a

-- | Function determining the validity of a state.
--
-- Useful when parameters are constrained and when calculating the prior or
-- likelihood functions is slow, but testing for validity is fast.
--
-- Also the validity of the state may depend on parameters untouched by the
-- Hamiltonian proposal.
type Validate a = a -> Bool

-- | Parameter mass matrix.
--
-- The masses roughly describe how reluctant the particles move through the
-- state space. If a parameter has higher mass, the momentum in this direction
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
--   to 0.0, and trust the tuning algorithm (see 'HTune') to find the correct
--   values.
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
-- NOTE: To avoid errors, the left bound has an additional hard minimum of 1,
-- and the right bound is required to be larger equal than the left bound.
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

-- | Parameters and functions required by the Hamilton Monte Carlo proposal.
data HSpec a = HSpec
  { -- | The sample state is used for error checks and to calculate the dimension
    -- of the proposal.
    hSample :: a,
    -- | Extract values to be manipulated by the Hamiltonian proposal from the
    -- state.
    hToVector :: a -> Positions,
    -- | Put those values back into the state.
    hFromVectorWith :: a -> Positions -> a,
    hGradient :: Gradient a,
    hMaybeValidate :: Maybe (Validate a)
  }

checkHSpecWith :: Eq a => HTuningSpec -> HSpec a -> Maybe String
checkHSpecWith tspec (HSpec x toVec fromVec _ _)
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

-- Internal. Symmetric, inverted mass matrix.
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

-- | Internal; Leapfrog integrator (also used by NUTS proposal).
leapfrog ::
  Gradient Positions ->
  Maybe (Validate Positions) ->
  HMassesInv ->
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  Positions ->
  Momenta ->
  -- Maybe (Positions', Momenta'); fail if state is not valid.
  Maybe (Positions, Momenta)
leapfrog grad mVal hMassesInv l eps theta phi = do
  let -- The first half step of the momenta.
      phiHalf = leapfrogStepMomenta (0.5 * eps) grad theta phi
  -- L-1 full steps. This gives the positions theta_{L-1}, and the momenta
  -- phi_{L-1/2}.
  (thetaLM1, phiLM1Half) <- go (l - 1) (Just (theta, phiHalf))
  -- The last full step of the positions.
  thetaL <- valF $ leapfrogStepPositions hMassesInv eps thetaLM1 phiLM1Half
  let -- The last half step of the momenta.
      --
      -- TODO (high): Since the gradient is evaluated at the final thetaL, one
      -- could use 'grad'', which also calculates the value of the posterior. Of
      -- course, this is only possible if the proposal data type is changed and
      -- one can provide a posterior (see comments in "Nuts").
      phiL = leapfrogStepMomenta (0.5 * eps) grad thetaL phiLM1Half
  return (thetaL, phiL)
  where
    valF x = case mVal of
      Nothing -> Just x
      Just f -> if f x then Just x else Nothing
    go _ Nothing = Nothing
    go n (Just (t, p))
      | n <= 0 = Just (t, p)
      | otherwise =
          let t' = leapfrogStepPositions hMassesInv eps t p
              p' = leapfrogStepMomenta eps grad t' p
              r = (,p') <$> valF t'
           in go (n - 1) r

leapfrogStepMomenta ::
  LeapfrogScalingFactor ->
  Gradient Positions ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New momenta.
  Momenta
leapfrogStepMomenta eps grad theta phi = phi + L.scale eps (grad theta)

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
  HTuningSpec ->
  HSpec a ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (ProposalSimple a)
hamiltonianSimpleWithTuningParameters tspec hspec t ts = do
  tspec' <- hTuningParametersToTuningSpec tspec t ts
  pure $ hamiltonianSimple tspec' hspec

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianSimple'.
hamiltonianSimpleWithMemoizedCovariance ::
  HTuningSpec ->
  HSpec a ->
  HData ->
  ProposalSimple a
hamiltonianSimpleWithMemoizedCovariance tspec hspec dt x g = do
  phi <- generateMomenta mu masses g
  lRan <- uniformR (lL, lR) g
  eRan <- uniformR (eL, eR) g
  case leapfrog gradientVec mValVec massesInv lRan eRan theta phi of
    -- TODO (high): Use the gradient to calculate the prior and the posterior.
    Nothing -> return (x, 0.0, 1.0, Nothing, Nothing)
    Just (theta', phi') ->
      let -- Prior of momenta.
          prPhi = logDensityMultivariateNormal mu massesInv logDetMasses phi
          prPhi' = logDensityMultivariateNormal mu massesInv logDetMasses phi'
          kernelR = prPhi' / prPhi
       in -- NOTE: For example, Neal page 12: In order for the Hamiltonian proposal
          -- to be in detailed balance, the momenta have to be negated before
          -- proposing the new value. That is, the negated momenta would guide the
          -- chain back to the previous state. However, we are only interested in
          -- the positions, and are not even storing the momenta.
          --
          -- TODO (high): Use the gradient to calculate the prior and the posterior.
          return (fromVec x theta', kernelR, 1.0, Nothing, Nothing)
  where
    (HTuningSpec masses l e _) = tspec
    (HSpec _ toVec fromVec gradient mVal) = hspec
    theta = toVec x
    lL = maximum [1 :: Int, floor $ (0.8 :: Double) * fromIntegral l]
    lR = maximum [lL, ceiling $ (1.2 :: Double) * fromIntegral l]
    eL = 0.8 * e
    eR = 1.2 * e
    (HData mu massesInv logDetMasses) = dt
    -- Vectorize the gradient and validation functions.
    gradientVec = toVec . gradient . fromVec x
    mValVec = mVal >>= (\f -> return $ f . fromVec x)

hamiltonianSimple :: HTuningSpec -> HSpec a -> ProposalSimple a
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
  Eq a =>
  HTuningSpec ->
  HSpec a ->
  PName ->
  PWeight ->
  Proposal a
hamiltonian tspec hspec n w = case checkHSpecWith tspec hspec of
  Just err -> error err
  Nothing ->
    let desc = PDescription "Hamiltonian Monte Carlo (HMC)"
        toVec = hToVector hspec
        dim = (L.size $ toVec $ hSample hspec)
        pDim = PSpecial dim 0.65
        ts = massesToTuningParameters (hMasses tspec)
        ps = hamiltonianSimple tspec hspec
        hamiltonianWith = Proposal n desc PSlow pDim w ps
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
            let tuner = Tuner 1.0 tFun ts tFunAux (hamiltonianSimpleWithTuningParameters tspec hspec)
             in hamiltonianWith $ Just tuner
