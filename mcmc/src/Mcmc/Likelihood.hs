-- |
-- Module      :  Mcmc.Likelihood
-- Description :  Types and convenience functions for computing likelihoods
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Wed Mar  3 11:39:04 2021.
module Mcmc.Likelihood
  ( Likelihood,
    LikelihoodFunction,
    LikelihoodFunctionG,
    noLikelihood,
  )
where

import Numeric.Log

-- | Likelihood values are stored in log domain.
type Likelihood = Log Double

-- | Likelihood function.
type LikelihoodFunction a = a -> Log Double

-- | Generalized likelihood function.
type LikelihoodFunctionG a b = a -> Log b

-- | Flat likelihood function. Useful for testing and debugging.
noLikelihood :: (RealFloat b) => LikelihoodFunctionG a b
noLikelihood = const 1.0
{-# SPECIALIZE noLikelihood :: LikelihoodFunction Double #-}
