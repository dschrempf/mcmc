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
--
-- This module provides prior distributions.
module Prior
  ( -- * Prior distributions
    uniformWith,
    normalWith,
    exponentialWith,
    gammaWith,
    branchesWith,

    -- * Auxiliary functions
    product',
  )
where

import ELynx.Data.Tree
import Mcmc
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Exponential
import Statistics.Distribution.Gamma
import Statistics.Distribution.Normal

-- | Uniform prior in [a, b].
uniformWith ::
  -- | Lower bound a.
  Double ->
  -- | Upper bound b.
  Double ->
  Double ->
  Log Double
uniformWith a b x
  | x <= a = pzero
  | x >= b = pzero
  | otherwise = Exp 0

-- | Normal distributed prior.
normalWith ::
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  Double ->
  Log Double
normalWith m s x = Exp $ logDensity d x
  where d = normalDistr m s


-- | Exponentail prior.
exponentialWith ::
  -- | Rate.
  Double ->
  Double ->
  Log Double
exponentialWith l x = Exp $ logDensity d x
  where
    d = exponential l

-- | Gamma distributed prior.
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

-- | Intelligent product that stops when encountering a zero.
product' :: [Log Double] -> Log Double
product' [] = 0
product' [x] = x
product' (x : xs)
  | x == 0 = 0
  | otherwise = x * product' xs
