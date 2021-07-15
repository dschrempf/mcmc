{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveGeneric #-}

-- |
-- Module      :  Statistics.Distribution.TruncatedNormal
-- Description :  Truncated normal distribution
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Sep  2 15:46:49 2020.
--
-- The truncated normal distribution is the probability distribution derived
-- from that of a normally distributed random variable by bounding the random
-- variable from either below or above (or both).
module Statistics.Distribution.TruncatedNormal
  ( TruncatedNormalDistribution,

    -- * Constructors
    truncatedNormalDistr,
  )
where

import Data.Data
import GHC.Generics
import Mcmc.Statistics.Types
import Numeric.MathFunctions.Constants
import Numeric.SpecFunctions
import qualified Statistics.Distribution as D

-- | The truncated normal distribution.
--
-- See https://en.wikipedia.org/wiki/Truncated_normal_distribution.
data TruncatedNormalDistribution = TND
  { mean :: {-# UNPACK #-} !(Mean Double),
    stdDev :: {-# UNPACK #-} !(StandardDeviation Double),
    lowerB :: {-# UNPACK #-} !(LowerBoundary Double),
    upperB :: {-# UNPACK #-} !(UpperBoundary Double),
    tndPhi2Alpha :: {-# UNPACK #-} !Double,
    tndDenom :: {-# UNPACK #-} !Double
  }
  deriving (Eq, Typeable, Data, Generic)

instance D.Distribution TruncatedNormalDistribution where
  cumulative = cumulative

instance D.ContDistr TruncatedNormalDistribution where
  density = density
  quantile = quantile

-- | Create a truncated normal distribution from parameters.
truncatedNormalDistr ::
  Mean Double ->
  StandardDeviation Double ->
  LowerBoundary Double ->
  UpperBoundary Double ->
  TruncatedNormalDistribution
truncatedNormalDistr m s a b
  | s <= 0 = error "truncatedNormalDistr: Standard deviation must be positive."
  | a >= b =
    error $
      "truncatedNormalDistr: Lower bound "
        <> show a
        <> " is equal or larger upper bound "
        <> show b
        <> "."
  | a > m || b < m =
    error $
      "truncatedNormalDistr: Mean "
        ++ show m
        ++ " is out of bounds ("
        ++ show a
        ++ ","
        ++ show b
        ++ ")."
  | otherwise = TND m s a b phi2Alpha (phi2 beta - phi2Alpha)
  where
    alpha = (a - m) / s
    beta = (b - m) / s
    phi2Alpha = phi2 alpha

phi1 :: Double -> Double
phi1 x = recip m_sqrt_2_pi * exp (- 0.5 * x * x)

phi2 :: Double -> Double
phi2 x = 0.5 * (1 + erf (x * m_1_sqrt_2))

getXi :: TruncatedNormalDistribution -> Double -> Double
getXi d x = (x - m) / s
  where
    m = mean d
    s = stdDev d

density :: TruncatedNormalDistribution -> Double -> Double
density d x
  | x < lowerB d = 0
  | x > upperB d = 0
  | otherwise = recip s * recip z * phi1 xi
  where
    s = stdDev d
    z = tndDenom d
    xi = getXi d x

cumulative :: TruncatedNormalDistribution -> Double -> Double
cumulative d x
  | x < lowerB d = 0
  | x > upperB d = 1
  | otherwise = recip z * (phi2 xi - phi2Alpha)
  where
    xi = getXi d x
    phi2Alpha = tndPhi2Alpha d
    z = tndDenom d

quantile :: TruncatedNormalDistribution -> Double -> Double
quantile d p
  | p == 0 = lowerB d
  | p == 1 = upperB d
  | p < 0 =
    error "quantile: p is smaller than 0."
  | p > 1 =
    error "quantile: p is larger than 0."
  | otherwise = invErf val * m_sqrt_2 * s + m
  where
    m = mean d
    s = stdDev d
    z = tndDenom d
    phi2Alpha = tndPhi2Alpha d
    val = 2 * (p * z + phi2Alpha) - 1
