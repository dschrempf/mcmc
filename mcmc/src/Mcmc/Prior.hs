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
module Mcmc.Prior
  ( Prior,
    PriorG,
    PriorFunction,
    PriorFunctionG,

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
import Mcmc.Internal.Gamma
import Mcmc.Statistics.Types
import Numeric.Log
import qualified Statistics.Distribution as S
import qualified Statistics.Distribution.Poisson as S

-- | Prior values are stored in log domain.
type Prior = PriorG Double

-- | Generalized prior.
type PriorG a = Log a

-- | Prior function.
type PriorFunction a = PriorFunctionG a Double

-- | Generalized prior function.
type PriorFunctionG a b = a -> PriorG b

-- | Flat prior function. Useful for testing and debugging.
noPrior :: RealFloat b => PriorFunctionG a b
noPrior = const 1.0
{-# SPECIALIZE noPrior :: PriorFunction Double #-}

-- | Improper uniform prior; strictly greater than a given value.
greaterThan :: RealFloat a => LowerBoundary a -> PriorFunctionG a a
greaterThan a x
  | x > a = 1.0
  | otherwise = 0.0
{-# SPECIALIZE greaterThan :: Double -> PriorFunction Double #-}

-- | Improper uniform prior; strictly greater than zero.
positive :: RealFloat a => PriorFunctionG a a
positive = greaterThan 0
{-# SPECIALIZE positive :: PriorFunction Double #-}

-- | Improper uniform prior; strictly less than a given value.
lessThan :: RealFloat a => UpperBoundary a -> PriorFunctionG a a
lessThan a x
  | x < a = 1.0
  | otherwise = 0.0
{-# SPECIALIZE lessThan :: Double -> PriorFunction Double #-}

-- | Improper uniform prior; strictly less than zero.
negative :: RealFloat a => PriorFunctionG a a
negative = lessThan 0
{-# SPECIALIZE negative :: PriorFunction Double #-}

-- | Exponential distributed prior.
exponential :: RealFloat a => Rate a -> PriorFunctionG a a
exponential l x = ll * Exp (negate l * x)
  where
    ll = Exp $ log l
{-# SPECIALIZE exponential :: Double -> PriorFunction Double #-}

-- | Gamma distributed prior.
gamma :: RealFloat a => Shape a -> Scale a -> PriorFunctionG a a
gamma k t x
  | x <= 0 = 0.0
  | otherwise = Exp $ log x * (k - 1) - (x / t) - logGammaG k - log t * k
{-# SPECIALIZE gamma :: Double -> Double -> PriorFunction Double #-}

-- | See 'gamma' but parametrized using mean and variance.
gammaMeanVariance :: RealFloat a => Mean a -> Variance a -> PriorFunctionG a a
gammaMeanVariance m v = gamma k t
  where
    (k, t) = gammaMeanVarianceToShapeScale m v
{-# SPECIALIZE gammaMeanVariance :: Double -> Double -> PriorFunction Double #-}

-- | Gamma disstributed prior with given shape and mean 1.0.
gammaMeanOne :: RealFloat a => Shape a -> PriorFunctionG a a
gammaMeanOne k = gamma k (recip k)
{-# SPECIALIZE gammaMeanOne :: Double -> PriorFunction Double #-}

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
gammaShapeScaleToMeanVariance :: Num a => Shape a -> Scale a -> (Mean a, Variance a)
gammaShapeScaleToMeanVariance k t = let m = k * t in (m, m * t)
{-# SPECIALIZE gammaShapeScaleToMeanVariance :: Double  -> Double -> (Double, Double) #-}

-- | Calculate shape and scale of the gamma distribution given the mean and
-- the variance.
gammaMeanVarianceToShapeScale :: Fractional a => Mean a -> Variance a -> (Shape a, Scale a)
gammaMeanVarianceToShapeScale m v = (m * m / v, v / m)
{-# SPECIALIZE gammaMeanVarianceToShapeScale :: Double  -> Double -> (Double, Double) #-}

mLnSqrt2Pi :: RealFloat a => a
mLnSqrt2Pi = 0.9189385332046727417803297364056176398613974736377834128171
{-# INLINE mLnSqrt2Pi #-}

-- | Normal distributed prior.
normal :: RealFloat a => Mean a -> StandardDeviation a -> PriorFunctionG a a
normal m s x = Exp $ (- xm * xm / (2 * s * s)) - denom
  where
    xm = x - m
    denom = mLnSqrt2Pi + log s
{-# SPECIALIZE normal :: Double -> Double -> PriorFunction Double #-}

-- | Uniform prior on [a, b].
uniform :: RealFloat a => LowerBoundary a -> UpperBoundary a -> PriorFunctionG a a
uniform a b x
  | x < a = 0.0
  | x > b = 0.0
  | otherwise = 1.0
{-# SPECIALIZE uniform :: Double -> Double -> PriorFunction Double #-}

-- | Poisson distributed prior.
poisson :: Rate Double -> PriorFunction Int
poisson l = Exp . S.logProbability d
  where
    d = S.poisson l

-- | Intelligent product that stops when encountering a zero.
--
-- Use with care because the elements are checked for positiveness, and this can
-- take some time if the list is long and does not contain any zeroes.
product' :: RealFloat a => [Log a] -> Log a
product' = fromMaybe 0 . prodM
{-# SPECIALIZE product' :: [Log Double] -> Log Double #-}

-- The type could be generalized to any MonadPlus Integer
prodM :: RealFloat a => [Log a] -> Maybe (Log a)
prodM = foldM (\ !acc x -> (acc * x) <$ guard (acc /= 0)) 1
{-# SPECIALIZE prodM :: [Log Double] -> Maybe (Log Double) #-}
