{-# LANGUAGE BangPatterns #-}

-- |
-- Module      :  Prior
-- Description :  Types and convenience functions for computing priors
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 13:26:14 2020.
--
-- Specialized prior functions. For the generalized versions, see
-- "Mcmc.Prior.General".
module Mcmc.Prior
  ( Prior,
    PriorFunction,

    -- * Improper priors
    noPrior,
    greaterThan,
    positive,
    lessThan,
    negative,

    -- * Continuous priors
    exponential,
    gamma,
    gammaMeanVariance,
    gammaMeanOne,
    gammaShapeScaleToMeanVariance,
    gammaMeanVarianceToShapeScale,
    normal,
    uniform,

    -- * Discrete priors
    poisson,

    -- * Auxiliary functions
    product',
  )
where

import Control.Monad
import Data.Maybe (fromMaybe)
import Mcmc.Statistics.Types
import Numeric.Log
import qualified Statistics.Distribution as S
import qualified Statistics.Distribution.Exponential as S
import qualified Statistics.Distribution.Gamma as S
import qualified Statistics.Distribution.Normal as S
import qualified Statistics.Distribution.Poisson as S

-- | Prior values are stored in log domain.
type Prior = Log Double

-- | Prior function.
type PriorFunction a = a -> Prior

-- | Flat prior function. Useful for testing and debugging.
noPrior :: PriorFunction a
noPrior = const 1.0

-- | Improper uniform prior; strictly greater than a given value.
greaterThan :: LowerBoundary -> PriorFunction Double
greaterThan a x
  | x <= a = 0
  | otherwise = 1

-- | Improper uniform prior; strictly greater than zero.
positive :: PriorFunction Double
positive = greaterThan 0

-- | Improper uniform prior; strictly less than a given value.
lessThan :: UpperBoundary -> PriorFunction Double
lessThan b x
  | x >= b = 0
  | otherwise = 1

-- | Improper uniform prior; strictly less than zero.
negative :: PriorFunction Double
negative = lessThan 0

-- | Exponential distributed prior.
exponential :: Rate -> PriorFunction Double
exponential l = Exp . S.logDensity d
  where
    d = S.exponential l

-- | Gamma distributed prior.
gamma :: Shape -> Scale -> PriorFunction Double
gamma k t = Exp . S.logDensity d
  where
    d = S.gammaDistr k t

-- | See 'gamma' but parametrized using mean and variance.
gammaMeanVariance :: Mean -> Variance -> PriorFunction Double
gammaMeanVariance m v = Exp . S.logDensity d
  where
    (k, th) = gammaMeanVarianceToShapeScale m v
    d = S.gammaDistr k th

-- | Gamma disstributed prior with given shape and mean 1.0.
gammaMeanOne :: Shape -> PriorFunction Double
gammaMeanOne k = Exp . S.logDensity d
  where
    d = S.gammaDistr k (recip k)

-- The mean and variance of the gamma distribution are
--
-- m = k*t
--
-- v = k*t*t
--
-- Hence, the shape and scale are
--
-- k = m^2/v
--
-- t = v/m

-- | Calculate mean and variance of the gamma distribution given the shape and
-- the scale.
gammaShapeScaleToMeanVariance :: Num a => ShapeG a -> ScaleG a -> (MeanG a, VarianceG a)
gammaShapeScaleToMeanVariance k t = let m = k * t in (m, m * t)

-- | Calculate shape and scale of the gamma distribution given the mean and
-- the variance.
gammaMeanVarianceToShapeScale :: Fractional a => MeanG a -> VarianceG a -> (ShapeG a, ScaleG a)
gammaMeanVarianceToShapeScale m v = (m * m / v, v / m)

-- | Normal distributed prior.
normal :: Mean -> StandardDeviation -> PriorFunction Double
normal m s = Exp . S.logDensity d
  where
    d = S.normalDistr m s

-- | Uniform prior on [a, b].
uniform :: LowerBoundary -> UpperBoundary -> PriorFunction Double
uniform a b x
  | x < a = 0
  | x > b = 0
  | otherwise = 1.0

-- | Poisson distributed prior.
poisson :: Rate -> PriorFunction Int
poisson l = Exp . S.logProbability d
  where
    d = S.poisson l

-- | Intelligent product that stops when encountering a zero.
--
-- Use with care because the elements are checked for positiveness, and this can
-- take some time if the list is long and does not contain any zeroes.
product' :: [Log Double] -> Log Double
product' = fromMaybe 0 . prodM

-- The type could be generalized to any MonadPlus Integer
prodM :: [Log Double] -> Maybe (Log Double)
prodM = foldM (\ !acc x -> (acc * x) <$ guard (acc /= 0)) 1
