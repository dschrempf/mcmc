{-# LANGUAGE TemplateHaskell #-}

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
    beta,
  )
where

import Data.Aeson
import Data.Aeson.TH
import qualified Data.Vector.Unboxed as V
import Mcmc.Proposal
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Beta
import Statistics.Distribution.Dirichlet

-- import Debug.Trace

-- | An element of a simplex.
--
-- A vector of non-negative values summing to one.
--
-- The nomenclature is not very consistent, because a K-dimensional simplex is
-- usually considered to be the set containing all @K@-dimensional vectors with
-- non-negative elements that sum to 1.0. However, I couldn't come up with a
-- better name. Maybe @SimplexElement@, but that was too long.
newtype Simplex = SimplexUnsafe {toVector :: V.Vector Double}
  deriving (Eq, Show)

$(deriveJSON defaultOptions ''Simplex)

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
  | otherwise = Right $ SimplexUnsafe v

-- | Create the uniform element of the K-dimensional simplex.
--
-- Set all values to \(1/D\).
simplexUniform :: Int -> Simplex
simplexUniform k = either error id $ simplexFromVector $ V.replicate k (1.0 / fromIntegral k)

-- Tuning function is inverted (high alpha means small steps).
getTuningFunction :: Double -> (Double -> Double)
getTuningFunction t = (/ t'')
  where
    -- Start with small steps.
    t' = t / 10000
    -- Extremely small tuning parameters lead to numeric overflow. The square
    -- root pulls the tuning parameter closer to 1.0. However, overflow may
    -- still occur (the involved Gamma functions grow faster than the
    -- exponential). I did not observe numeric underflow in my tests.
    t'' = sqrt t'

-- The tuning parameter is proportional to the inverted mean of the shape
-- parameter values.
--
-- The values determining the proposal size have been set using an example
-- analysis. They are good values for this analysis, but may fail for other
-- analyses.
dirichletSimple :: Double -> ProposalSimple Simplex
dirichletSimple t (SimplexUnsafe xs) g = do
  -- If @t@ is high and above 1.0, the parameter vector will be low, and the
  -- variance will be high. If @t@ is low and below 1.0, the parameter vector
  -- will be high, and the Dirichlet distribution will be very concentrated with
  -- low variance.
  let ddXs = either error id $ dirichletDistribution $ V.map tf xs
  -- traceShowM $ V.map tf xs
  ys <- dirichletSample ddXs g
  -- traceShowM ys
  -- Have to check if parameters are valid (because zeroes do occur).
  let eitherDdYs = dirichletDistribution $ V.map tf ys
  let mhRatio = case eitherDdYs of
        -- Set ratio to 0; so that the proposal will not be accepted.
        Left _ -> 0
        Right ddYs -> dirichletDensity ddYs xs / dirichletDensity ddXs ys
  -- traceShowM mhRatio
  return (SimplexUnsafe ys, mhRatio)
  where
    tf = getTuningFunction t

-- | Dirichlet proposal on a simplex.
--
-- For a given element of a K-dimensional simplex, propose a new element of the
-- K-dimensional simplex. The new element is sampled from the multivariate
-- Dirichlet distribution with parameter vector being the old element of the
-- simplex.
--
-- The tuning parameter is used to determine the concentration of the Dirichlet
-- distribution: the lower the tuning parameter, the higher the concentration.
--
-- This proposal may have low acceptance ratios. In this case, please see the
-- coordinate wise 'beta' proposal.
dirichlet :: String -> Int -> Bool -> Proposal Simplex
dirichlet = createProposal "Dirichlet" dirichletSimple

-- The tuning parameter is the inverted mean of the shape values.
--
-- The values determining the proposal size have been set using an example
-- analysis. They are good values for this analysis, but may fail for other
-- analyses.
--
-- See also the 'dirichlet' proposal.
betaSimple :: Int -> Double -> ProposalSimple Simplex
betaSimple i t (SimplexUnsafe xs) g = do
  -- Shape parameters of beta distribution. Do not assume that the sum of the
  -- elements of 'xs' is 1.0, because then repeated proposals let the sum of the
  -- vector diverge.
  let aX = xI
      bX = xsSum - xI
      bdXI = betaDistr (tf aX) (tf bX)
  yI <- genContVar bdXI g
  -- Shape parameters of beta distribution.
  let aY = yI
      bY = 1.0 - yI
      eitherBdYI = betaDistrE (tf aY) (tf bY)
  -- See 'dirichlet', which has the same construct.
  let mhRatio = case eitherBdYI of
        Nothing -> 0
        Just bdYI -> Exp $ logDensity bdYI xI - logDensity bdXI yI
  -- Construct new vector.
  let
    nf x = x * bY / bX
    ys = V.generate (V.length xs) (\j -> if i==j then yI else nf (xs V.! j))
  return (either error id $ simplexFromVector ys, mhRatio)
  where
    xI = xs V.! i
    xsSum = V.sum xs
    tf = getTuningFunction t

-- | Beta proposal on a specific coordinate @i@ on a simplex.
--
-- For a given element of a K-dimensional simplex, propose a new element of the
-- K-dimensional simplex. The coordinate @i@ of the new element is sampled from
-- the beta distribution. The other coordinates are normalized such that the
-- values sum to 1.0. The parameters of the beta distribution are chosen such
-- that the expected value of the beta distribution is the value of the old
-- coordinate.
--
-- The tuning parameter is used to determine the concentration of the beta
-- distribution: the lower the tuning parameter, the higher the concentration.
--
-- See also the 'dirichlet' proposal.
--
-- No "out of bounds" checks are performed during compile time. Run time errors
-- can occur if @i@ is negative, or if @i-1@ is larger than the length of the
-- element vector of the simplex.
beta :: Int -> String -> Int -> Bool -> Proposal Simplex
beta i = createProposal description (betaSimple i)
  where description = "Beta; coordinate: " ++ show i