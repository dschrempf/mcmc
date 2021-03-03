-- |
-- Module      :  Mcmc.Likelihood
-- Description :  Types and convenience functions for computing likelihoods
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Wed Mar  3 11:39:04 2021.
module Mcmc.Likelihood
  ( LikelihoodFunction,
    noLikelihood,
  )
where

import Numeric.Log

-- | Likelihood function.
type LikelihoodFunction a = a -> Log Double

-- | Flat likelihood function. Useful for testing and debugging.
noLikelihood :: LikelihoodFunction a
noLikelihood = const 1.0
