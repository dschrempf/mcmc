-- |
-- Module      :  Mcmc.Statistics.Types
-- Description :  Types indicating properties of distributions
-- Copyright   :  2021 Dominik Schrempf
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
    Rate,
    Dimension,
    Size,
    LowerBoundary,
    UpperBoundary,
  )
where

-- | Mean of a distribution.
type Mean a = a

-- | Standard deviation of a distribution.
type StandardDeviation a = a

-- | Variance of a distribution.
type Variance a = a

-- | Shape of a distribution.
type Shape a = a

-- | Scale of a distribution.
type Scale a = a

-- | Rate of a distribution.
type Rate a = a

-- | Dimension of a distribution.
type Dimension = Int

-- | Size of a distribution.
--
-- For example, the size of the interval of the uniform distribution.
type Size = Double

-- | Lower boundary of a distribution.
type LowerBoundary a = a

-- | Upper boundary of a distribution.
type UpperBoundary a = a
