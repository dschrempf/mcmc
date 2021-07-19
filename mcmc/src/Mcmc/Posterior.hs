-- |
-- Module      :  Mcmc.Posterior
-- Description :  Types for posterior values and functions
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri May 28 12:26:35 2021.
module Mcmc.Posterior
  ( Posterior,
    PosteriorG,
    PosteriorFunction,
    PosteriorFunctionG,
  )
where

import Numeric.Log

-- | Posterior values are stored in log domain.
type Posterior = PosteriorG Double

-- | Generalized posterior.
type PosteriorG a = Log a

-- | Posterior function.
type PosteriorFunction a = PosteriorFunctionG a Double

-- | Generalized posterior function.
type PosteriorFunctionG a b = a -> PosteriorG b
