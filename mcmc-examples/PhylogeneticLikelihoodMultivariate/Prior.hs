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
    uniform,
    normal,
    exponential,
    gamma,
    branchesWith,

    -- * Auxiliary functions
    product',
  )
where

import ELynx.Data.Tree
import Numeric.Log
import qualified Statistics.Distribution as S
import qualified Statistics.Distribution.Exponential as S
import qualified Statistics.Distribution.Gamma as S
import qualified Statistics.Distribution.Normal as S

-- | Larger than 0.
positive :: Double -> Log Double
positive x | x <= 0 = 0
           | otherwise = 1

-- | Uniform prior in [a, b].
uniform ::
  -- | Lower bound a.
  Double ->
  -- | Upper bound b.
  Double ->
  Double ->
  Log Double
uniform a b x
  | x <= a = 0
  | x >= b = 0
  | otherwise = Exp 0

-- | Normal distributed prior.
normal ::
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  Double ->
  Log Double
normal m s x = Exp $ S.logDensity d x
  where d = S.normalDistr m s

-- | Exponentail prior.
exponential ::
  -- | Rate.
  Double ->
  Double ->
  Log Double
exponential l x = Exp $ S.logDensity d x
  where
    d = S.exponential l

-- | Gamma distributed prior.
gamma ::
  -- | Shape.
  Double ->
  -- | Scale.
  Double ->
  Double ->
  Log Double
gamma k t x = Exp $ S.logDensity d x
  where
    d = S.gammaDistr k t

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
