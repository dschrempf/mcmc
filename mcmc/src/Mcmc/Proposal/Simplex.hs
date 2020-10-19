-- |
-- Module      :  Mcmc.Proposal.Simplex
-- Description :  Proposals on simplices
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Oct 19 15:32:31 2020.
module Mcmc.Proposal.Simplex
  ( -- * Simplices
    Simplex (toVector),
    simplexUniform,
    simplexFromVector,
    simplexFromVectorNormalize,
  )
where

import Control.Monad.Primitive
import qualified Data.Vector.Unboxed as V
import Numeric.Log hiding (sum)
import Numeric.SpecFunctions
import System.Random.MWC (Gen)
import System.Random.MWC.Distributions (gamma)

newtype Simplex = Simplex {toVector :: V.Vector Double}
  deriving (Eq, Show)

-- | Set all values to \(1/n\).
simplexUniform :: Int -> Simplex
simplexUniform n | n > 0 = Simplex $ V.replicate n (1.0 / fromIntegral n)
                 | otherwise = error "simplexUniform: Integer must be strictly positive."

-- Tolerance.
eps :: Double
eps = 1e-14

-- Check if vector is normalized with tolerance 'eps'.
isNormalized :: V.Vector Double -> Bool
isNormalized v
  | abs (V.sum v - 1.0) > eps = False
  | otherwise = True

-- | Create a simplex from a vector.
--
-- - Call 'error' if vector is not normalized.
-- - Call 'error' if vector is empty.
simplexFromVector :: V.Vector Double -> Simplex
simplexFromVector v
  | V.null v = error "simplexFromVector: Vector is empty."
  | not (isNormalized v) = error "simplexFromVector: Vector is not normalized."
  | otherwise = Simplex v

-- | Create a simplex from a vector. The vector is normalized to sum to 1.0.
--
-- Call 'error' if vector is empty.
simplexFromVectorNormalize :: V.Vector Double -> Simplex
simplexFromVectorNormalize v
  | V.null v = error "simplexFromVectorNormalize: Vector is empty."
  | otherwise = Simplex $ V.map (/ s) v
  where
    s = V.sum v

-- Dirichlet density for a given parameter and value vector.
dirichletDensity :: V.Vector Double -> V.Vector Double -> Log Double
dirichletDensity as xs
  | nAs /= nXs = error "dirichletDensity: Parameter and value vectors have different length."
  | nAs == 0 = error "dirichletDensity: Parameter vector is empty."
  | nXs == 0 = error "dirichletDensity: Value vector is empty."
  | otherwise =
    if not (isNormalized xs)
      then 0
      else Exp $ logDenominator - logNominator + logXsPow
  where
    nAs = V.length as
    nXs = V.length xs
    logNominator = V.sum $ V.map logGamma as
    logDenominator = logGamma (V.sum as)
    logXsPow = V.sum $ V.zipWith (\a x -> log $ x ** (a - 1.0)) as xs

-- Sample a value vector from the Dirichlet distribution with given parameter
-- vector and a generator.
dirichletSample :: PrimMonad m => V.Vector Double -> Gen (PrimState m) -> m (V.Vector Double)
dirichletSample as g = do
  ys <- V.mapM (\a -> gamma a 1.0 g) as
  let s = V.sum ys
  return $ V.map (/ s) ys
