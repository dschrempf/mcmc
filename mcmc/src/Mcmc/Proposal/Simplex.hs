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
  ( -- * Elements of simplices
    Simplex (toVector),
    simplexUniform,
    simplexFromVector,

    -- * Proposals on simplices
    dirichlet,
  )
where

-- TODO: SimplexElementScale (?).

import qualified Data.Vector.Unboxed as V
import Mcmc.Proposal
import Statistics.Distribution.Dirichlet

-- | An element of a simplex.
--
-- A vector of non-negative values summing to one.
--
-- The nomenclature is not very consistent, because a K-dimensional simplex is
-- usually considered to be the set containing all @K@-dimensional vectors with
-- non-negative elements that sum to 1.0. However, I couldn't come up with a
-- better name. Maybe @SimplexElement@, but that was too long.
newtype Simplex = Simplex {toVector :: V.Vector Double}
  deriving (Eq, Show)

-- Tolerance.
eps :: Double
eps = 1e-14

-- Check if vector is normalized with tolerance 'eps'.
isNormalized :: V.Vector Double -> Bool
isNormalized v
  | abs (V.sum v - 1.0) > eps = False
  | otherwise = True

-- Check if vector contains negative elements.
isNegative :: V.Vector Double -> Bool
isNegative = V.any (< 0)

-- | Ensure that the value vector is an element of a simplex.
--
-- Return 'Left' if:
-- - The value vector is empty.
-- - The value vector contains negative elements.
-- - The value vector is not normalized.
simplexFromVector :: V.Vector Double -> Either String Simplex
simplexFromVector v
  | V.null v = Left "simplexFromVector: Vector is empty."
  | isNegative v = Left "simplexFromVector: Vector contains negative elements."
  | not (isNormalized v) = Left "simplexFromVector: Vector is not normalized."
  | otherwise = Right $ Simplex v

-- | Create the uniform element of the K-dimensional simplex.
--
-- Set all values to \(1/D\).
simplexUniform :: Int -> Simplex
simplexUniform k = either error id $ simplexFromVector $ V.replicate k (1.0 / fromIntegral k)

-- The tuning parameter is the inverted mean of all alpha values.
dirichletSimple :: Double -> ProposalSimple Simplex
dirichletSimple t (Simplex xs) g = do
  -- If @t@ is high and above 1.0, the parameter vector will be low, and the
  -- variance will be high. If @t@ is low and below 1.0, the parameter vector
  -- will be high, and the Dirichlet distribution will be very concentrated with
  -- low variance.
  let ddXs = either error id $ dirichletDistribution $ V.map (/t) xs
  ys <- dirichletSample ddXs g
  let ddYs = either error id $ dirichletDistribution $ V.map (/t) ys
      densityXY = dirichletDensity ddXs ys
      densityYX = dirichletDensity ddYs xs
  return (either error id $ simplexFromVector ys, densityYX / densityXY)

-- | Dirichlet proposal.
--
-- For a given element of a K-dimensional simplex, propose a new element of the
-- K-dimensional simplex. The new element is sampled from the multivariate
-- Dirichlet distribution with parameter vector being the old element of the
-- simplex. The tuning parameter is used to determine the concentration of the
-- Dirichlet distribution: the lower the tuning parameter, the higher the
-- concentration.
dirichlet :: String -> Int -> Bool -> Proposal Simplex
dirichlet = createProposal dirichletSimple
