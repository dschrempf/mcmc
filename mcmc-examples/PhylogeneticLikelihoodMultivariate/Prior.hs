-- |
-- Module      :  Prior
-- Description :  Functions to compute standard priors
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 13:26:14 2020.
module Prior
  ( exponentialWith,
    gammaWith,
    branchesWith,
  )
where

import Data.Foldable
import ELynx.Data.Tree
-- import Mcmc
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Exponential
import Statistics.Distribution.Gamma

-- | Prior with exponential distribution.
exponentialWith ::
  -- | Rate.
  Double ->
  Double ->
  Log Double
exponentialWith l x = Exp $ logDensity d x
  where
    d = exponential l

-- | Prior with gamma distribution.
gammaWith ::
  -- | Shape.
  Double ->
  -- | Scale.
  Double ->
  Double ->
  Log Double
gammaWith k t x = Exp $ logDensity d x
  where
    d = gammaDistr k t

-- | Branch length prior with given distribution.
--
-- Root branch is ignored!
branchesWith ::
  -- | Branch prior distribution.
  (Double -> Log Double) ->
  Tree Double a ->
  Log Double
branchesWith f = product . map f . tail . branches
