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
    positive,
    uniformWith,
    normalWith,
    exponentialWith,
    gammaWith,
    birthAndDeathWith,
    branchesWith,

    -- * Auxiliary functions
    product',
  )
where

import ELynx.Data.Tree
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Exponential
import Statistics.Distribution.Gamma
import Statistics.Distribution.Normal

-- | Larger than 0.
positive :: Double -> Log Double
positive x | x <= 0 = 0
           | otherwise = 1

-- | Uniform prior in [a, b].
uniformWith ::
  -- | Lower bound a.
  Double ->
  -- | Upper bound b.
  Double ->
  Double ->
  Log Double
uniformWith a b x
  | x <= a = 0
  | x >= b = 0
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

-- TODO: First the heights have to be accessible on the tree. Then use the birth
-- and death distribution from elynx-tree with origin or root height.
--
-- Question: Is the birth death distribution of the point process adequate, if
-- the topology is given?
--
-- How do we handle the critical, or nearly critical process? Is it fast enough
-- to combine the three distributions into one, and check during calculation of
-- the density? Or is it better to check for the three different possibilities
-- here, and then use the adequate distribution (status quo)?

-- | Birth and death prior.
birthAndDeathWith ::
  -- | Birth rate lambda.
  Double ->
  -- | Death rate mu.
  Double ->
  Tree e a ->
  Log Double
birthAndDeathWith l m t = undefined

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
