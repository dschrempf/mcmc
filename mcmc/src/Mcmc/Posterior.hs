-- |
-- Module      :  Mcmc.Posterior
-- Description :  Types for posterior values and functions
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri May 28 12:26:35 2021.
module Mcmc.Posterior
  ( Posterior,
    PosteriorFunction,
    PosteriorFunctionG,
  )
where

import Numeric.Log

-- | Posterior values are stored in log domain.
type Posterior = Log Double

-- | Posterior function.
type PosteriorFunction a = a -> Log Double

-- | Generalized posterior function.
type PosteriorFunctionG a b = a -> Log b
