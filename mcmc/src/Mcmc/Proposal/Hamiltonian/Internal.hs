-- |
-- Module      :  Mcmc.Proposal.Hamiltonian.Internal
-- Description :  Internal definitions related to Hamiltonian dynamics
-- Copyright   :  2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Thu Jun  9 15:12:39 2022.
module Mcmc.Proposal.Hamiltonian.Internal
  ( -- * Tuning
    findReasonableEpsilon,
    massesToTuningParameters,
    tuningParametersToMasses,
    tuneDiagonalMassesOnly,
    tuneAllMasses,

    -- * Structure of state
    checkHStructureWith,

    -- * Hamiltonian dynamics
    Mu,
    MassesInv,
    HData (..),
    getHData,
    generateMomenta,
    exponentialKineticEnergy,

    -- * Leapfrog integrator
    Target,
    leapfrog,
  )
where

import qualified Data.Vector as VB
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Mcmc.Proposal
import Mcmc.Proposal.Hamiltonian.Common
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import qualified Statistics.Covariance as S
import qualified Statistics.Function as S
import qualified Statistics.Sample as S
import System.Random.MWC
import System.Random.Stateful

-- See Algorithm 4 in Matthew D. Hoffman, Andrew Gelman (2014) The No-U-Turn
-- Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo, Journal
-- of Machine Learning Research.
findReasonableEpsilon ::
  StatefulGen g m =>
  Target ->
  Masses ->
  HData ->
  Positions ->
  g ->
  m LeapfrogScalingFactor
findReasonableEpsilon t ms (HData mu msInv) q g = do
  p <- generateMomenta mu ms g
  case leapfrog t msInv 1 eI q p of
    Nothing -> undefined
    Just (_, p', prQ, prQ') -> do
      let expEKin = exponentialKineticEnergy msInv p
          expEKin' = exponentialKineticEnergy msInv p'
          rI :: Double
          rI = exp $ ln $ prQ' * expEKin' / (prQ * expEKin)
          a :: Double
          a = if rI > 0.5 then 1 else (-1)
          go e r =
            if r ** a > 2 ** (negate a)
              then case leapfrog t msInv 1 e q p of
                Nothing -> e
                Just (_, p'', _, prQ'') ->
                  let expEKin'' = exponentialKineticEnergy msInv p''
                      r' :: Double
                      r' = exp $ ln $ prQ'' * expEKin'' / (prQ * expEKin)
                      e' = (2 ** a) * e
                   in go e' r'
              else e
      pure $ go eI rI
  where
    eI = 1.0

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

-- | Internal.
massesToTuningParameters :: Masses -> AuxiliaryTuningParameters
massesToTuningParameters = VB.convert . L.flatten . L.unSym

-- | Internal.
tuningParametersToMasses ::
  -- Dimension of the mass matrix.
  Int ->
  AuxiliaryTuningParameters ->
  Masses
tuningParametersToMasses d = L.trustSym . L.reshape d . VB.convert

-- | Internal.
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
-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneDiagonalMassesOnly dim toVec xs ts
  -- If not enough data is available, do not tune.
  | VB.length xs < samplesMinDiagonal = ts
  | otherwise =
      -- Replace the diagonal.
      massesToTuningParameters $
        L.trustSym $
          massesOld - L.diag massesDiagonalOld + L.diag massesDiagonalNew
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

-- | Internal.
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
-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
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
    sigmaNormalized = L.unSym $ either error fst $ S.graphicalLasso 0.1 xsNormalized
    sigma = S.rescaleSWith ss sigmaNormalized
    massesNew = L.inv sigma

-- | Internal.
checkHStructureWith :: Eq (s Double) => Masses -> HStructure s -> Maybe String
checkHStructureWith ms (HStructure x toVec fromVec)
  | fromVec x xVec /= x = eWith "'fromVectorWith x (toVector x) /= x' for sample state."
  | L.size xVec /= nrows = eWith "Mass matrix and 'toVector x' have different sizes for sample state."
  | otherwise = Nothing
  where
    eWith m = Just $ "checkHStructureWith: " <> m
    nrows = L.rows $ L.unSym ms
    xVec = toVec x

-- | Internal. Mean vector containing zeroes. We save this vector because it is
-- required when sampling from the multivariate normal distribution.
type Mu = L.Vector Double

-- | Internal. Symmetric, inverted mass matrix.
type MassesInv = L.Herm Double

-- | Internal data type containing memoized values.
data HData = HData
  { hMu :: Mu,
    hMassesInv :: MassesInv
  }
  deriving (Show)

-- | Internal. Compute inverted mass matrix.
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

-- | Internal. Generate momenta for a new iteration.
generateMomenta ::
  StatefulGen g m =>
  Mu ->
  Masses ->
  g ->
  m Momenta
generateMomenta mu masses gen = do
  seed <- uniformM gen
  let momenta = L.gaussianSample seed 1 mu masses
  return $ L.flatten momenta

-- TODO (medium): Use a sparse matrix approach for the log density of the
-- multivariate normal, similar to McmcDate.

-- | Internal. Compute exponent of kinetic energy.
exponentialKineticEnergy ::
  MassesInv ->
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

-- | Internal; Leapfrog integrator.
leapfrog ::
  Target ->
  MassesInv ->
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
leapfrog tF msInv l eps theta phi = do
  let -- The first half step of the momenta.
      (x, phiHalf) = leapfrogStepMomenta (0.5 * eps) tF theta phi
  -- L-1 full steps for positions and momenta. This gives the positions
  -- theta_{L-1}, and the momenta phi_{L-1/2}.
  (thetaLM1, phiLM1Half) <- go (l - 1) $ Just $ (theta, phiHalf)
  -- The last full step of the positions.
  let thetaL = leapfrogStepPositions msInv eps thetaLM1 phiLM1Half
  -- The last half step of the momenta.
  let (x', phiL) = leapfrogStepMomenta (0.5 * eps) tF thetaL phiLM1Half
  return (thetaL, phiL, x, x')
  where
    go _ Nothing = Nothing
    go n (Just (t, p))
      | n <= 0 = Just (t, p)
      | otherwise =
          let t' = leapfrogStepPositions msInv eps t p
              (x, p') = leapfrogStepMomenta eps tF t' p
           in if x > 0.0
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
  MassesInv ->
  LeapfrogScalingFactor ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New positions.
  Positions
leapfrogStepPositions msInv eps theta phi = theta + (L.unSym msInvEps L.#> phi)
  where
    msInvEps = L.scale eps msInv
