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
-- - The Hamiltonian proposal acts on a vector of storable 'Values'. Functions
--   converting the state to and from this vector have to be provided. See
--   'HSettings'.
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
  ( Values,
    Gradient,
    Validate,
    Masses,
    LeapfrogTrajectoryLength,
    LeapfrogScalingFactor,
    HTuneLeapfrog (..),
    HTuneMasses (..),
    HTune (..),
    HSettings (..),
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

-- TODO: No-U-turn sampler (NUTS). Ameliorates necessity to determine the
-- leapfrog trajectory length L. (I think this is a necessary extension.)

-- TODO: Riemannian adaptation: State-dependent mass matrix. (Seems a little bit
-- of an overkill.)

-- | The Hamiltonian proposal acts on a vector of floating point values.
type Values = L.Vector Double

-- | Gradient of the log posterior function.
--
-- The gradient has to be provided for the complete state. The reason is that
-- the gradient may change if parameters untouched by the Hamiltonian proposal
-- are altered by other proposals.
type Gradient a = a -> a

-- | Function validating the state.
--
-- Useful if parameters are constrained.
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

-- Internal. Values; target state containing parameters.
type Positions = Values

-- Internal. Momenta of the parameters.
type Momenta = L.Vector Double

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

-- | Tuning settings.
data HTune = HTune HTuneLeapfrog HTuneMasses
  deriving (Eq, Show)

-- | Specifications of the Hamilton Monte Carlo proposal.
data HSettings a = HSettings
  { -- | Extract values to be manipulated by the Hamiltonian proposal from the
    -- state.
    hToVector :: a -> Values,
    -- | Put those values back into the state.
    hFromVectorWith :: a -> Values -> a,
    hGradient :: Gradient a,
    hMaybeValidate :: Maybe (Validate a),
    hMasses :: Masses,
    hLeapfrogTrajectoryLength :: LeapfrogTrajectoryLength,
    hLeapfrogScalingFactor :: LeapfrogScalingFactor,
    hTune :: HTune
  }

checkHSettings :: Eq a => a -> HSettings a -> Maybe String
checkHSettings x (HSettings toVec fromVec _ _ masses l eps _)
  | any (<= 0) diagonalMasses = Just "checkHSettings: Some diagonal entries of the mass matrix are zero or negative."
  | nrows /= ncols = Just "checkHSettings: Mass matrix is not square."
  | fromVec x xVec /= x = Just "checkHSettings: 'fromVectorWith x (toVector x) /= x' for sample state."
  | L.size xVec /= nrows = Just "checkHSettings: Mass matrix has different size than 'toVector x', where x is sample state."
  | l < 1 = Just "checkHSettings: Leapfrog trajectory length is zero or negative."
  | eps <= 0 = Just "checkHSettings: Leapfrog scaling factor is zero or negative."
  | otherwise = Nothing
  where
    ms = L.unSym masses
    diagonalMasses = L.toList $ L.takeDiag ms
    nrows = L.rows ms
    ncols = L.cols ms
    xVec = toVec x

-- Internal. Mean vector containing zeroes.
type HMu = L.Vector Double

-- Internal. Symmetric, inverted mass matrix.
type HMassesInv = L.Herm Double

-- Internal. Symmetric, inverted mass matrix scaled with the leapfrog step size
-- epsilon.
type HMassesInvEps = L.Herm Double

-- Internal. Logarithm of the determinant of the mass matrix.
type HLogDetMasses = Double

-- Internal data type containing memoized values.
data HData = HData
  { _hMu :: HMu,
    _hMassesInv :: HMassesInv,
    _hMassesInvEps :: HMassesInvEps,
    _hLogDetMasses :: HLogDetMasses
  }

-- Call 'error' if the determinant of the covariance matrix is negative.
getHData :: HSettings a -> HData
getHData s =
  -- The multivariate normal distribution requires a positive definite matrix
  -- with positive determinant.
  if sign == 1.0
    then HData mu massesInvH massesInvEpsH logDetMasses
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
    eps = hLeapfrogScalingFactor s
    massesInvEpsH = L.scale eps massesInvH

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

leapfrog ::
  Gradient Positions ->
  Maybe (Validate Positions) ->
  HMassesInvEps ->
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  Positions ->
  Momenta ->
  -- Maybe (Positions', Momenta'); fail if state is not valid.
  Maybe (Positions, Momenta)
leapfrog grad mVal hMassesInvEps l eps theta phi = do
  let -- The first half step of the momenta.
      phiHalf = leapfrogStepMomenta (0.5 * eps) grad theta phi
  -- L-1 full steps. This gives the positions theta_{L-1}, and the momenta
  -- phi_{L-1/2}.
  (thetaLM1, phiLM1Half) <- go (l - 1) (Just (theta, phiHalf))
  -- The last full step of the positions.
  thetaL <- valF $ leapfrogStepPositions hMassesInvEps thetaLM1 phiLM1Half
  let -- The last half step of the momenta.
      phiL = leapfrogStepMomenta (0.5 * eps) grad thetaL phiLM1Half
  return (thetaL, phiL)
  where
    valF x = case mVal of
      Nothing -> Just x
      Just f -> if f x then Just x else Nothing
    go _ Nothing = Nothing
    go 0 (Just (t, p)) = Just (t, p)
    go n (Just (t, p)) =
      let t' = leapfrogStepPositions hMassesInvEps t p
          p' = leapfrogStepMomenta eps grad t' p
          r = (,p') <$> valF t'
       in go (n - 1) r

leapfrogStepMomenta ::
  -- Step size (usually a multiple of the leapfrog scaling factor).
  Double ->
  Gradient Positions ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New momenta.
  Momenta
leapfrogStepMomenta eps grad theta phi = phi + L.scale eps (grad theta)

leapfrogStepPositions ::
  HMassesInvEps ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New positions.
  Positions
leapfrogStepPositions hMassesInvEps theta phi = theta + (L.unSym hMassesInvEps L.#> phi)

massesToTuningParameters :: Masses -> AuxiliaryTuningParameters
massesToTuningParameters = VB.convert . L.flatten . L.unSym

tuningParametersToMasses ::
  -- Dimension of the mass matrix.
  Int ->
  AuxiliaryTuningParameters ->
  Masses
tuningParametersToMasses d = L.trustSym . L.reshape d . VB.convert

hTuningParametersToSettings ::
  HSettings a ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (HSettings a)
hTuningParametersToSettings s t ts
  | nTsNotOK =
      Left "hTuningParametersToSettings: Auxiliary variables do not have correct dimension."
  | otherwise =
      Right $
        s
          { hMasses = msTuned,
            hLeapfrogTrajectoryLength = lTuned,
            hLeapfrogScalingFactor = eTuned
          }
  where
    ms = hMasses s
    d = L.rows $ L.unSym ms
    l = hLeapfrogTrajectoryLength s
    e = hLeapfrogScalingFactor s
    (HTune tlf tms) = hTune s
    nTsNotOK =
      let nTs = VU.length ts
       in case tms of
            HNoTuneMasses -> nTs /= 0
            _ -> nTs /= d * d
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
  HSettings a ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (ProposalSimple a)
hamiltonianSimpleWithTuningParameters s t ts =
  hamiltonianSimple <$> hTuningParametersToSettings s t ts

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianSimple'.
hamiltonianSimpleWithMemoizedCovariance ::
  HSettings a ->
  HData ->
  ProposalSimple a
hamiltonianSimpleWithMemoizedCovariance st dt x g = do
  phi <- generateMomenta mu masses g
  lRan <- uniformR (lL, lR) g
  eRan <- uniformR (eL, eR) g
  case leapfrog gradientVec mValVec massesInvEps lRan eRan theta phi of
    Nothing -> return (x, 0.0, 1.0)
    Just (theta', phi') ->
      let -- Prior of momenta.
          prPhi = logDensityMultivariateNormal mu massesInv logDetMasses phi
          -- NOTE: Neal page 12: In order for the proposal to be in detailed
          -- balance, the momenta have to be negated before proposing the new
          -- value. This is not required here since the prior involves a
          -- multivariate normal distribution with means 0.
          prPhi' = logDensityMultivariateNormal mu massesInv logDetMasses phi'
          kernelR = prPhi' / prPhi
       in return (fromVec x theta', kernelR, 1.0)
  where
    (HSettings toVec fromVec gradient mVal masses l e _) = st
    theta = toVec x
    lL = maximum [1 :: Int, floor $ (0.8 :: Double) * fromIntegral l]
    lR = maximum [lL, ceiling $ (1.2 :: Double) * fromIntegral l]
    eL = 0.8 * e
    eR = 1.2 * e
    (HData mu massesInv massesInvEps logDetMasses) = dt
    -- Vectorize the gradient and validation functions.
    gradientVec = toVec . gradient . fromVec x
    mValVec = mVal >>= (\f -> return $ f . fromVec x)

hamiltonianSimple ::
  HSettings a ->
  ProposalSimple a
hamiltonianSimple s = hamiltonianSimpleWithMemoizedCovariance s hd
  where
    hd = getHData s

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
  -- | The sample state is used for error checks and to calculate the dimension
  -- of the proposal.
  a ->
  HSettings a ->
  PName ->
  PWeight ->
  Proposal a
hamiltonian x s n w = case checkHSettings x s of
  Just err -> error err
  Nothing ->
    let desc = PDescription "Hamiltonian Monte Carlo (HMC)"
        toVec = hToVector s
        dim = (L.size $ toVec x)
        pDim = PSpecial dim 0.65
        ts = massesToTuningParameters (hMasses s)
        ps = hamiltonianSimple s
        hamiltonianWith = Proposal n desc pDim w ps
        tSet@(HTune tlf tms) = hTune s
        tFun = case tlf of
          HNoTuneLeapfrog -> noTuningFunction
          HTuneLeapfrog -> defaultTuningFunctionWith pDim
        tFunAux = case tms of
          HNoTuneMasses -> noAuxiliaryTuningFunction
          HTuneDiagonalMassesOnly -> tuneDiagonalMassesOnly dim toVec
          HTuneAllMasses -> tuneAllMasses dim toVec
     in case tSet of
          (HTune HNoTuneLeapfrog HNoTuneMasses) -> hamiltonianWith Nothing
          _ ->
            let tuner = Tuner 1.0 tFun ts tFunAux (hamiltonianSimpleWithTuningParameters s)
             in hamiltonianWith $ Just tuner
