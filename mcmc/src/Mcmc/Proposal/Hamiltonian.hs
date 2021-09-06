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
-- - The implementation assumes the existence of the gradient. Like so, the user
--   can use automatic or manual differentiation, depending on the problem at
--   hand.
--
-- - The state needs to be list like or 'Traversable' so that the structure of
--   the state space is available. A 'Traversable' constraint on the data type
--   is nice because it is more general than, for example, a list, and
--   user-defined data structures can be used.
--
-- - The state needs to have a zip-like 'Applicative' instance so that
-- - matrix/vector operations can be performed.
module Mcmc.Proposal.Hamiltonian
  ( Gradient,
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

import Data.Foldable
-- TODO: Remove this!
import qualified Data.Matrix as M
import Data.Maybe
import Data.Traversable
import qualified Data.Vector as VB
import Mcmc.Prior
import Mcmc.Proposal
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import Numeric.MathFunctions.Constants
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified Statistics.Function as S
import qualified Statistics.Sample as S
import System.Random.MWC

-- TODO: At the moment, the HMC proposal is agnostic of the prior and
-- likelihood, that is, the posterior function. This means, that it cannot know
-- when it reaches a point with zero posterior probability. This also affects
-- restricted or constrained parameters. See Gelman p. 303.
--
-- Mon Sep 6 02:22:38 PM CEST 2021: This is not entirely true anymore, see
-- 'Validate'.

-- TODO: No-U-turn sampler.

-- TODO: Riemannian adaptation.

-- | Gradient of the log posterior function.
type Gradient = L.Vector Double -> L.Vector Double

-- | Function validating the state.
--
-- Useful if parameters are constrained.
type Validate = L.Vector Double -> Bool

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
--   to 0.0, and trust the tuning algorithm (see 'HTuneMassesAndLeapfrog') to
--   find the correct values.
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
-- NOTE: To avoid errors, the left bound has an additional hard minimum of 1,
-- and the right bound is required to be larger equal than the left bound.
--
-- Usually set to 10, but larger values may be desirable.
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

-- Target state containing parameters.
type Positions = L.Vector Double

-- Momenta of the parameters.
type Momenta = L.Vector Double

-- | Tuning of leapfrog parameters:
--
-- We expect that the larger the leapfrog step size the larger the proposal step
-- size and the lower the acceptance ratio. Consequently, if the acceptance rate
-- is too low, the leapfrog step size is decreased and vice versa. Further, the
-- leapfrog trajectory length is scaled such that the product of the leapfrog
-- step size and trajectory length stays constant.
data HTuneLeapfrog = HNoTuneLeapfrog | HTuneLeapfrog
  deriving (Eq, Show)

-- | Tuning of masses:
--
-- The variances of all parameters of the posterior distribution obtained over
-- the last auto tuning interval is calculated and the masses are amended using
-- the old masses and the inverted variances. If, for a specific coordinate, the
-- sample size is too low, or if the calculated variance is out of predefined
-- bounds, the mass of the affected position is not changed.
data HTuneMasses = HNoTuneMasses | HTuneDiagonalMassesOnly | HTuneAllMasses
  deriving (Eq, Show)

-- | Tuning settings.
data HTune = HTune HTuneLeapfrog HTuneMasses
  deriving (Eq, Show)

-- | Specifications for Hamilton Monte Carlo proposal.
data HSettings a = HSettings
  { hToVector :: a -> L.Vector Double,
    hFromVector :: L.Vector Double -> a,
    hGradient :: Gradient,
    hMaybeValidate :: Maybe Validate,
    hMasses :: Masses,
    hLeapfrogTrajectoryLength :: LeapfrogTrajectoryLength,
    hLeapfrogScalingFactor :: LeapfrogScalingFactor,
    hTune :: HTune
  }

checkHSettings :: a -> HSettings a -> Either String ()
checkHSettings x (HSettings toV fromV _ _ masses l eps _)
  | any (<= 0) diagonals = Left "checkHSettings: Some diagonal entries of the mass matrix are zero or negative."
  | nrows /= ncols = Left "checkHSettings: Mass matrix is not square."
  | fromV xVec /= x = Left "checkHSettings: 'fromVector . toVector' is not 'id' for sample state."
  | L.size xVec /= nrows = Left "checkHSettings: Mass matrix has different size than 'toVector x', where x is sample state."
  | l < 1 = Left "checkHSettings: Leapfrog trajectory length is zero or negative."
  | eps <= 0 = Left "checkHSettings: Leapfrog scaling factor is zero or negative."
  | otherwise = Right ()
  where
    ms = L.unSym masses
    diagonals = L.toList $ L.takeDiag ms
    nrows = L.rows ms
    ncols = L.cols ms
    xVec = toV x

generateMomenta ::
  -- Mean vector (zeroes). Provided so that it does not have to be created.
  L.Vector Double ->
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
  -- Log of determinant of covariance matrix.
  Double ->
  -- Value vector.
  L.Vector Double ->
  Log Double
logDensityMultivariateNormal mu sigmaInvH logSigmaDet xs =
  Exp $ c + (-0.5) * (logSigmaDet + ((dxs L.<# sigmaInv) L.<.> dxs))
  where
    dxs = xs - mu
    k = fromIntegral $ L.size mu
    c = negate $ m_ln_sqrt_2_pi * k
    sigmaInv = L.unSym sigmaInvH

leapfrog ::
  Gradient ->
  Maybe Validate ->
  Masses ->
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  Positions ->
  Momenta ->
  -- Maybe (Positions', Momenta').
  Maybe (Positions, Momenta)
leapfrog grad mVal masses l eps theta phi = do
  let -- The first half step of the momenta.
      phiHalf = leapfrogStepMomenta 0.5 eps grad theta phi
  -- L-1 full steps. This gives the positions theta_{L-1}, and the momenta
  -- phi_{L-1/2}.
  (thetaLM1, phiLM1Half) <- go (l - 1) (Just (theta, phiHalf))
  -- The last full step of the positions.
  thetaL <- valF $ leapfrogStepPositions eps masses thetaLM1 phiLM1Half
  let -- The last half step of the momenta.
      phiL = leapfrogStepMomenta 0.5 eps grad thetaL phiLM1Half
  return (thetaL, phiL)
  where
    valF x = case mVal of
      Nothing -> Just x
      Just f -> if f x then Just x else Nothing
    go _ Nothing = Nothing
    go 0 (Just (t, p)) = Just (t, p)
    go n (Just (t, p)) =
      let t' = leapfrogStepPositions eps masses t p
          p' = leapfrogStepMomenta 1.0 eps grad t' p
          r = (,p') <$> valF t'
       in go (n - 1) r

leapfrogStepMomenta ::
  -- Size of step (half or full step).
  Double ->
  LeapfrogScalingFactor ->
  Gradient ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New momenta.
  Momenta
leapfrogStepMomenta xi eps grad theta phi = phi <+. ((xi * eps) .* grad theta)
  where
    (<+.) :: Applicative f => f (Maybe Double) -> f Double -> f (Maybe Double)
    (<+.) xs ys = f <$> xs <*> ys
    f Nothing _ = Nothing
    f (Just x) y = Just $ x + y

leapfrogStepPositions ::
  LeapfrogScalingFactor ->
  Masses ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  Positions
-- The arguments are flipped to encounter the maybe momentum.
leapfrogStepPositions eps masses theta phi = theta <+. (mScaledReversed .*> phi)
  where
    (<+.) :: Applicative f => f Double -> f (Maybe Double) -> f Double
    (<+.) xs ys = f <$> xs <*> ys
    f x Nothing = x
    f x (Just y) = x + y
    mScaledReversed = (fmap . fmap) ((* eps) . (** (-1))) masses
    (.*>) :: Applicative f => f (Maybe Double) -> f (Maybe Double) -> f (Maybe Double)
    (.*>) xs ys = g <$> xs <*> ys
    g (Just x) (Just y) = Just $ x * y
    g Nothing Nothing = Nothing
    g _ _ = error "leapfrogStepPositions: Got just a mass and no momentum or the other way around."

-- Scalar-vector multiplication.
(.*) :: Applicative f => Double -> f Double -> f Double
(.*) x ys = (* x) <$> ys

-- NOTE: Fixed parameters without mass have a tuning parameter of NaN.
massesToTuningParameters :: Masses -> AuxiliaryTuningParameters
massesToTuningParameters = VB.fromList . map (fromMaybe nan) . toList
  where
    nan = 0 / 0

-- We need the structure in order to fill it with the given parameters.
tuningParametersToMasses ::
  AuxiliaryTuningParameters ->
  Masses ->
  Either String Masses
tuningParametersToMasses xs ms =
  if null xs'
    then sequenceA msE
    else Left "tuningParametersToMasses: Too many values."
  where
    (xs', msE) = mapAccumL setValue (VB.toList xs) ms
    setValue [] _ = ([], Left "tuningParametersToMasses: Too few values.")
    -- NOTE: Recover fixed parameters and unset their mass.
    setValue (y : ys) _ = let y' = if isNaN y then Nothing else Just y in (ys, Right y')

hTuningParametersToSettings ::
  TuningParameter ->
  AuxiliaryTuningParameters ->
  HSettings a ->
  Either String (HSettings a)
hTuningParametersToSettings t ts (HSettings _ _ g v m l e tn) =
  if tn == HTuneMassesAndLeapfrog
    then case tuningParametersToMasses ts m of
      Left err -> Left err
      Right m' -> Right $ HSettings g v m' lTuned eTuned tn
    else Right $ HSettings g v m lTuned eTuned tn
  where
    -- The larger epsilon, the larger the proposal step size and the lower the
    -- expected acceptance ratio.
    --
    -- Further, we roughly keep \( L * \epsilon = 1.0 \). The equation is not
    -- correct, because we pull L closer to the original value to keep the
    -- runtime somewhat acceptable.
    lTuned = ceiling $ fromIntegral l / (t ** 0.9) :: Int
    eTuned = t * e

hamiltonianSimpleWithTuningParameters ::
  HSettings a ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (ProposalSimple a)
hamiltonianSimpleWithTuningParameters s t ts = case hTuningParametersToSettings t ts s of
  Left err -> Left err
  Right s' -> Right $ hamiltonianSimple s'

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianSimple'.
hamiltonianSimpleMemoizedCovariance ::
  HSettings a ->
  -- Mean vector (zeroes).
  L.Vector Double ->
  -- Inverted covariance matrix.
  L.Herm Double ->
  -- Log of determinant of covariance matrix.
  Double ->
  ProposalSimple a
hamiltonianSimpleMemoizedCovariance s mu sigmaInv logSigmaDet theta g = do
  phi <- generateMomenta mu masses g
  lRan <- uniformR (lL, lR) g
  eRan <- uniformR (eL, eR) g
  case leapfrog gradient mVal masses lRan eRan theta phi of
    Nothing -> return (theta, 0.0, 1.0)
    Just (theta', phi') ->
      let -- Prior of momenta.
          prPhi = logDensityMultivariateNormal mu sigmaInv logSigmaDet phi
          -- NOTE: Neal page 12: In order for the proposal to be in detailed
          -- balance, the momenta have to be negated before proposing the new
          -- value. This is not required here since the prior involves a
          -- multivariate normal distribution with means 0.
          prPhi' = logDensityMultivariateNormal mu sigmaInv logSigmaDet phi'
          kernelR = prPhi' / prPhi
       in return (theta', kernelR, 1.0)
  where
    (HSettings _ _ gradient mVal masses l e _) = s
    lL = maximum [1 :: Int, floor $ (0.8 :: Double) * fromIntegral l]
    lR = maximum [lL, ceiling $ (1.2 :: Double) * fromIntegral l]
    eL = 0.8 * e
    eR = 1.2 * e

hamiltonianSimple ::
  HSettings a ->
  ProposalSimple a
hamiltonianSimple s =
  if sign == 1.0
    then hamiltonianSimpleMemoizedCovariance s mu sigmaInvH logSigmaDet
    else error "hamiltonianSimple: Determinant of covariance matrix is negative?"
  where
    ms = hMasses s
    nrows = L.rows $ L.unSym ms
    mu = L.fromList $ replicate nrows 0.0
    (sigmaInv, (logSigmaDet, sign)) = L.invlndet $ L.unSym ms
    -- In theory we can trust that the matrix is symmetric here, because the
    -- inverse of a symmetric matrix is symmetric. However, one may want to
    -- implement a check anyways.
    sigmaInvH = L.trustSym sigmaInv

minVariance :: Double
minVariance = 1e-6

maxVariance :: Double
maxVariance = 1e6

minSamples :: Int
minSamples = 60

computeAuxiliaryTuningParameters ::
  -- TODO: This is now a boxed vector of storable vectors. Please check if a
  -- change is required.
  VB.Vector Positions ->
  AuxiliaryTuningParameters ->
  AuxiliaryTuningParameters
computeAuxiliaryTuningParameters xss ts =
  VB.zipWith (\t -> rescueWith t . calcSamplesAndVariance) ts xssT
  where
    -- TODO: Improve matrix transposition. USE HMATRIX.
    xssT = VB.fromList $ M.toColumns $ M.fromLists $ VB.toList $ VB.map toList xss
    calcSamplesAndVariance xs = (VB.length $ VB.uniq $ S.gsort xs, S.variance xs)
    rescueWith t (sampleSize, var) =
      if var < minVariance || maxVariance < var || sampleSize < minSamples
        then -- then traceShow ("Rescue with " <> show t) t
          t
        else
          let t' = sqrt (t * recip var)
           in -- in traceShow ("Old mass " <> show t <> " new mass " <> show t') t'
              t'

-- | Hamiltonian Monte Carlo proposal.
--
-- The 'Applicative' and 'Traversable' instances are used for element-wise
-- operations.
--
-- Assume a zip-like 'Applicative' instance so that cardinality remains
-- constant.
--
-- NOTE: The desired acceptance rate is 0.65, although the dimension of the
-- proposal is high.
--
-- NOTE: The speed of this proposal can change drastically when tuned because
-- the leapfrog trajectory length is changed.
hamiltonian ::
  -- | The sample state is used to calculate the dimension of the proposal.
  a ->
  HSettings a ->
  PName ->
  PWeight ->
  Proposal a
hamiltonian x s n w = case checkHSettings x s of
  Just err -> error err
  Nothing ->
    let desc = PDescription "Hamiltonian Monte Carlo (HMC)"
        dim = PSpecial (length x) 0.65
        ts = massesToTuningParameters (hMasses s)
        ps = hamiltonianSimple s
        p' = Proposal n desc dim w ps
        fT = defaultTuningFunction dim
        tS = hTune s
        fTs =
          if tS == HTuneMassesAndLeapfrog
            then computeAuxiliaryTuningParameters
            else \_ xs -> xs
     in case tS of
          (HTune HNoTuneLeapfrog HNoTuneMasses) -> p' Nothing
          _ -> p' $ Just $ Tuner 1.0 fT ts fTs (hamiltonianSimpleWithTuningParameters s)
