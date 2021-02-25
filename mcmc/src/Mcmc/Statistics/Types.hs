-- |
-- Module      :  Mcmc.Statistics.Types
-- Description :  Type synonyms
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Tue Feb 23 14:51:06 2021.
module Mcmc.Statistics.Types
  ( Mean,
    StandardDeviation,
    Variance,
    Shape,
    Scale,
    Dimension,
    Size,
    LowerBoundary,
    UpperBoundary,
  )
where

-- | Type synonym indicating the mean of a distribution.
type Mean = Double

-- | Type synonym indicating the standard deviation of a distribution.
type StandardDeviation = Double

-- | Type synonym indicating the variance of a distribution.
type Variance = Double

-- | Type synonym indicating the shape of a distribution.
type Shape = Double

-- | Type synonym indicating the scale of a distribution.
type Scale = Double

-- | Type synonym indicating the dimension of a distribution.
type Dimension = Int

-- | Type synonym indicating the size of a distribution.
--
-- For example, the interval of the uniform distribution.
type Size = Double

-- | Type synonym indicating the lower boundary.
type LowerBoundary = Double

-- | Type synonym indicating the upper boundary.
type UpperBoundary = Double
