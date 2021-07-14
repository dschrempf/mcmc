-- |
-- Module      :  Mcmc.Statistics.Types
-- Description :  Types indicating properties of distributions
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
    MeanG,
    StandardDeviation,
    StandardDeviationG,
    Variance,
    VarianceG,
    Shape,
    ShapeG,
    Scale,
    ScaleG,
    Rate,
    RateG,
    Dimension,
    Size,
    LowerBoundary,
    LowerBoundaryG,
    UpperBoundary,
    UpperBoundaryG,
  )
where

-- | Mean of a distribution.
type Mean = Double

-- | Generalized 'Mean'.
type MeanG a = a

-- | Standard deviation of a distribution.
type StandardDeviation = Double

-- | Generalized 'StandardDeviation'.
type StandardDeviationG a = a

-- | Variance of a distribution.
type Variance = Double

-- | Generalized 'Variance'.
type VarianceG a = a

-- | Shape of a distribution.
type Shape = Double

-- | Generalized 'Shape'.
type ShapeG a = a

-- | Scale of a distribution.
type Scale = Double

-- | Generalized 'Scale'
type ScaleG a = a

-- | Rate of a distribution.
type Rate = Double

-- | Generalized 'Rate'.
type RateG a = a

-- | Dimension of a distribution.
type Dimension = Int

-- | Size of a distribution.
--
-- For example, the size of the interval of the uniform distribution.
type Size = Double

-- | Lower boundary of a distribution.
type LowerBoundary = Double

-- | Generalized 'LowerBoundary'.
type LowerBoundaryG a = a

-- | Upper boundary of a distribution.
type UpperBoundary = Double

-- | Generalized 'UpperBoundary'.
type UpperBoundaryG a = a
