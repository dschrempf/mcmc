-- |
-- Module      :  Mcmc.Proposal.Hamiltonian.Masses
-- Description :  Mass matrices
-- Copyright   :  2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Tue Jun 14 10:09:24 2022.
module Mcmc.Proposal.Hamiltonian.Masses
  ( Mu,
    MassesI,
    toGMatrix,
    cleanMatrix,
    getMassesI,
    getMus,
    Dimension,
    vectorToMasses,
    massesToVector,

    -- * Tuning
    SmoothingParameter (..),
    tuneDiagonalMassesOnly,
    tuneAllMasses,
  )
where

import Data.Maybe
import qualified Data.Vector as VB
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Mcmc.Proposal.Hamiltonian.Common
import qualified Numeric.LinearAlgebra as L
import Numeric.Natural
import qualified Statistics.Covariance as S
import qualified Statistics.Function as S
import qualified Statistics.Sample as S

-- Mean vector containing zeroes. We save this vector because it is required
-- when sampling from the multivariate normal distribution.
type Mu = L.Vector Double

-- General, symmetric, inverted mass matrix.
type MassesI = L.GMatrix

-- Purge masses and inverted masses (excluding the diagonal) strictly smaller
-- than the precision.
--
-- If changed, also change help text of 'HTuneMasses', which is indirectly
-- affected via 'massMin', and 'massMax'.
precision :: Double
precision = 1e-8

isDiag :: L.Matrix Double -> Bool
isDiag xs = abs (sumDiag - sumFull) < precision
  where
    xsAbs = L.cmap abs xs
    sumDiag = L.sumElements (L.takeDiag xsAbs)
    sumFull = L.sumElements xsAbs

-- Consider a matrix sparse if less than (5 * number of rows) elements are
-- non-zero.
isSparse :: L.Matrix Double -> Bool
isSparse xs = nNonZero < fromIntegral nMax
  where
    n = L.rows xs
    m = min 5 n
    nMax = n * m
    f x = if abs x >= precision then 1 else 0 :: Double
    xsInd = L.cmap f xs
    nNonZero = L.sumElements xsInd

toAssocMatrix :: L.Matrix Double -> L.AssocMatrix
toAssocMatrix xs
  | n /= m = error "toAssocMatrix: Matrix not square."
  | otherwise =
      [ ((i, j), e)
        | i <- [0 .. (n - 1)],
          j <- [0 .. (n - 1)],
          let e = xs `L.atIndex` (i, j),
          abs e >= precision
      ]
  where
    n = L.rows xs
    m = L.cols xs

toGMatrix :: L.Matrix Double -> L.GMatrix
toGMatrix xs
  | n == 0 || m == 0 = error "toGMatrix: Matrix empty."
  | n /= m = error "toGMatrix: Matrix not square."
  | isDiag xs = L.mkDiagR n m $ L.takeDiag xs
  | isSparse xs = L.mkSparse $ toAssocMatrix xs
  | otherwise = L.mkDense xs
  where
    n = L.rows xs
    m = L.cols xs

-- Diagonal:
-- - NaN values are set to 1.
-- - Negative values are set to 1.
-- - Small elements are set to 'precision'.
--
-- Off-diagonal:
-- - NaN values are set to 0.
-- - Elements with absolute values strictly smaller than 'precision' are purged.
--
-- We are permissive with negative and NaN values because adequate masses are
-- crucial. The Hamiltonian algorithms also work when the masses are off.
cleanMatrix :: L.Matrix Double -> L.Matrix Double
cleanMatrix xs =
  L.diag (L.cmap cleanDiag xsDiag) + L.cmap cleanOffDiag xsOffDiag
  where
    xsDiag = L.takeDiag xs
    cleanDiag x
      | isNaN x = 1
      | x < 0 = 1
      -- The strict comparison is important.
      | x < precision = precision
      | otherwise = x
    xsOffDiag = xs - L.diag xsDiag
    cleanOffDiag x
      | isNaN x = 0
      -- The strict comparison is important.
      | abs x < precision = 0
      | otherwise = x

getMassesI :: L.Herm Double -> L.GMatrix
getMassesI xs
  | n == 0 || m == 0 = error "getMassesI: Matrix empty."
  | n /= m = error "getMassesI: Matrix not square."
  | sign /= 1.0 = error "getMassesI: Determinant of matrix is negative."
  | otherwise = toGMatrix $ cleanMatrix xsI
  where
    xs' = L.unSym xs
    n = L.rows xs'
    m = L.cols xs'
    (xsI, (_, sign)) = L.invlndet xs'

getMus :: Masses -> L.Vector Double
getMus xs
  | n == 0 || m == 0 = error "getMu: Matrix empty."
  | n /= m = error "getMu: Matrix not square."
  | otherwise = L.fromList $ replicate n 0.0
  where
    xs' = L.unSym xs
    n = L.rows xs'
    m = L.cols xs'

-- Dimension of the proposal.
type Dimension = Int

massesToVector :: Masses -> VU.Vector Double
massesToVector = VU.convert . L.flatten . L.unSym

vectorToMasses :: Dimension -> VU.Vector Double -> Masses
vectorToMasses d = L.trustSym . L.reshape d . VU.convert

-- If changed, also change help text of 'HTuneMasses'.
massMin :: Double
massMin = precision

-- If changed, also change help text of 'HTuneMasses'.
massMax :: Double
massMax = recip precision

-- Minimal number of unique samples required for tuning the diagonal entries of
-- the mass matrix.
--
-- NOTE: If changed, also change help text of 'HTuneMasses'.
samplesDiagonalMin :: Int
samplesDiagonalMin = 61

-- Minimal number of samples required for tuning all entries of the mass matrix.
--
-- NOTE: If changed, also change help text of 'HTuneMasses'.
samplesAllMinWith :: Dimension -> Int
samplesAllMinWith d = samplesDiagonalMin + max samplesDiagonalMin d

getSampleSize :: VS.Vector Double -> Int
getSampleSize = VS.length . VS.uniq . S.gsort

rescueAfter :: Double -> Double
rescueAfter y
  | massMin > abs y = signum y * massMin
  | abs y > massMax = signum y * massMax
  | otherwise = y

rescueBefore :: (Double -> Double -> Double) -> Double -> Double -> Double
rescueBefore f old new
  -- Be permissive with NaN and infinite values.
  | isNaN new = old
  | isInfinite new = old
  | otherwise = f old new

rescue :: (Double -> Double -> Double) -> Double -> Double -> Double
rescue f old new = rescueAfter $ rescueBefore f old new

-- -- I do not know where I got this function from, but it works pretty well!
-- interpolate :: Double -> Double -> Double
-- interpolate old new = interSqrt ** 2
--   where
--     interSqrt = recip 3 * (sqrt old + 2 * sqrt new)

-- | This parameter plays the same role as @m@ in the dual averaging algorithm.
-- In the beginning, when @m@ is zero, the newly estimated masses will take full
-- precedence over the current masses. After some auto tuning steps, when @m@ is
-- larger, the newly estimated masses influence the current masses only
-- slightly.
newtype SmoothingParameter = SmoothingParameter Natural

-- Another interpolation function I came up with. It is pretty cool, because
-- (similar to above) it interpolates the square roots, which is what we want.
interpolate' :: SmoothingParameter -> Double -> Double -> Double
interpolate' (SmoothingParameter m) oldSquared newSquared = finSign * (fin ** 2)
  where
    oldSign = signum oldSquared
    newSign = signum newSquared
    sqrt' = sqrt . abs
    old = oldSign * sqrt' oldSquared
    new = newSign * sqrt' newSquared
    -- The new mass will be the second last boundary. That is, the larger the
    -- number of bins is, the more informative the new mass is compared to the
    -- old mass.
    nbins = max (100 - fromIntegral m) 2
    fin = recip nbins * (old + (nbins - 1) * new)
    finSign = signum fin

getNewMassDiagonalWithSampleSize :: SmoothingParameter -> Int -> Double -> Double -> Double
getNewMassDiagonalWithSampleSize m sampleSize massOld massEstimate
  | sampleSize < samplesDiagonalMin = massOld
  -- Diagonal masses are variances which are strictly positive.
  | massEstimate <= 0 = massOld
  | otherwise = rescue (interpolate' m) massOld massEstimate

getNewMassDiagonal :: SmoothingParameter -> Double -> Double -> Double
getNewMassDiagonal m massOld massEstimate
  -- Diagonal masses are variances which are strictly positive.
  | massEstimate <= 0 = massOld
  | otherwise = rescue (interpolate' m) massOld massEstimate

getNewMassOffDiagonal :: SmoothingParameter -> Double -> Double -> Double
getNewMassOffDiagonal m = rescue (interpolate' m)

-- The Cholesky decomposition, which is performed when sampling new momenta with
-- 'generateMomenta', requires a positive definite covariance matrix. The
-- Graphical Lasso algorithm finds positive definite covariance matrices, but
-- sometimes positive definiteness is violated because of numerical errors.
-- Further, when non-diagonal masses are already non-zero, the tuning of
-- diagonal masses only may violate positive definiteness.
--
-- Find the closest positive definite matrix of a given matrix.
--
-- See https://gist.github.com/fasiha/fdb5cec2054e6f1c6ae35476045a0bbd.
findClosestPositiveDefiniteMatrix :: L.Matrix Double -> L.Matrix Double
findClosestPositiveDefiniteMatrix a
  | n == 0 || m == 0 = error "findClosestPositiveDefiniteMatrix: Matrix empty."
  | n /= m = error "findClosestPositiveDefiniteMatrix: Matrix not square."
  | isPositiveDefinite a = a
  | otherwise = go a3 1
  where
    n = L.rows a
    m = L.cols a
    b = L.unSym $ L.sym a
    (_, s, v) = L.svd b
    h = L.tr v L.<> (L.diag s L.<> v)
    a2 = L.scale 0.5 (b + h)
    a3 = L.unSym $ L.sym a2
    isPositiveDefinite = isJust . L.mbChol . L.trustSym
    --
    i = L.ident n
    -- See https://hackage.haskell.org/package/ieee754-0.8.0/docs/src/Numeric-IEEE.html#line-177.
    eps = 2.2204460492503131e-16
    go x k
      | isPositiveDefinite x = x
      | otherwise =
          let minEig = L.minElement $ L.cmap L.realPart $ L.eigenvalues x
              nu = negate minEig * (k ** 2) + eps
              x' = x + L.scale nu i
           in go x' (k + 1)

tuneDiagonalMassesOnly ::
  SmoothingParameter ->
  -- Conversion from value to vector.
  (a -> Positions) ->
  -- Value vector.
  VB.Vector a ->
  -- Old mass matrix, and inverted mass matrix.
  (Masses, MassesI) ->
  -- new mass matrix, and inverted mass matrix.
  (Masses, MassesI)
-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneDiagonalMassesOnly m toVec xs (ms, msI)
  -- If not enough data is available, do not tune.
  | VB.length xs < samplesDiagonalMin = (ms, msI)
  | dimState /= dimMs = error "tuneDiagonalMassesOnly: Dimension mismatch."
  -- Replace the diagonal.
  | otherwise =
      let msDirty = msOld - L.diag msDiagonalOld + L.diag msDiagonalNew
          -- Positive definite matrices are symmetric.
          ms' = L.trustSym $ findClosestPositiveDefiniteMatrix $ cleanMatrix msDirty
          msI' = getMassesI ms'
       in (ms', msI')
  where
    -- xs: Each element contains all parameters of one iteration.
    -- xs': Each element is a vector containing all parameters changed by the
    -- proposal of one iteration.
    xs' = VB.map toVec xs
    -- xs'': Matrix with each row containing all parameter values changed by the
    -- proposal of one iteration.
    xs'' = L.fromRows $ VB.toList xs'
    -- We can safely use 'VB.head' here since the length of 'xs' must be larger
    -- than 'samplesDiagonalMin'.
    dimState = VS.length $ VB.head xs'
    sampleSizes = VS.fromList $ map getSampleSize $ L.toColumns xs''
    msOld = L.unSym ms
    dimMs = L.rows msOld
    msDiagonalOld = L.takeDiag msOld
    msDiagonalEstimate = VS.fromList $ map (recip . S.variance) $ L.toColumns xs''
    msDiagonalNew =
      VS.zipWith3
        (getNewMassDiagonalWithSampleSize m)
        sampleSizes
        msDiagonalOld
        msDiagonalEstimate

-- This value was carefully tuned using the example "hamiltonian".
defaultGraphicalLassoPenalty :: Double
defaultGraphicalLassoPenalty = 0.3

interpolateElementWise :: SmoothingParameter -> L.Matrix Double -> L.Matrix Double -> L.Matrix Double
interpolateElementWise m old new
  | mO /= mN = err "different number of rows"
  | nO /= nN = err "different number of columns"
  | mO /= nO = err "not square"
  | mO < 1 = err "empty matrix"
  | otherwise = L.build (mO, nO) f
  where
    mO = L.rows old
    nO = L.cols old
    mN = L.rows new
    nN = L.cols new
    err msg = error $ "interpolateElementWise: " <> msg
    f iD jD =
      let -- This sucks a bit, because we need a function (e -> e -> e), and
          -- since the return type is Double, the indices are also Doubble.
          i = round iD
          j = round jD
          g = if i == j then getNewMassDiagonal m else getNewMassOffDiagonal m
          eO = old `L.atIndex` (i, j)
          eN = new `L.atIndex` (i, j)
       in g eO eN

tuneAllMasses ::
  SmoothingParameter ->
  -- Conversion from value to vector.
  (a -> Positions) ->
  -- Value vector.
  VB.Vector a ->
  -- Old mass matrix, and inverted mass matrix.
  (Masses, MassesI) ->
  -- New mass matrix, and inverted mass matrix.
  (Masses, MassesI)
-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneAllMasses m toVec xs (ms, msI)
  -- If not enough data is available, do not tune.
  | VB.length xs < samplesDiagonalMin = (ms, msI)
  -- If not enough data is available, only the diagonal masses are tuned.
  | VB.length xs < samplesAllMinWith dimMs = fallbackDiagonal
  | L.rank xs'' /= dimState = fallbackDiagonal
  | dimState /= dimMs = error "tuneAllMasses: Dimension mismatch."
  | otherwise = (msNew, msINew)
  where
    fallbackDiagonal = tuneDiagonalMassesOnly m toVec xs (ms, msI)
    -- xs: Each element contains all parameters of one iteration.
    -- xs': Each element is a vector containing all parameters changed by the
    -- proposal of one iteration.
    xs' = VB.map toVec xs
    -- xs'': Matrix with each row containing all parameter values changed by the
    -- proposal of one iteration.
    xs'' = L.fromRows $ VB.toList xs'
    -- We can safely use 'VB.head' here since the length of 'xs' must be larger
    -- than 'samplesDiagonalMin'.
    dimState = VS.length $ VB.head xs'
    dimMs = L.rows $ L.unSym ms
    (_, ss, xsNormalized) = S.scale xs''
    -- The first value is the covariance matrix sigma, which the inverted mass
    -- matrix (precision matrix). However, we interpolate the new mass matrix
    -- using the old one and the new estimate, so we have to recalculate the
    -- covariance matrix anyways.
    (_, precNormalized) =
      either error id $
        S.graphicalLasso defaultGraphicalLassoPenalty xsNormalized
    ms' = S.rescalePWith ss (L.unSym precNormalized)
    -- Clean NaNs, infinities; ensure positive definiteness. The masses should
    -- be positive definite, but sometimes they happen to be not because of
    -- numerical errors.
    msNewDirty = interpolateElementWise m (L.unSym ms) ms'
    -- Positive definite matrices are symmetric.
    msNew = L.trustSym $ findClosestPositiveDefiniteMatrix $ cleanMatrix msNewDirty
    msINew = getMassesI msNew
