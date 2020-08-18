-- |
-- Module      :  Probability
-- Description :  Convenience functions to compute priors or posteriors
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 13:26:14 2020.
--
-- Convenience functions to compute priors or posteriors.
module Mcmc.Probability
  ( -- * Continuous probability density functions
    positive,
    uniform,
    normal,
    exponential,
    gamma,

    -- * Discrete probability mass functions
    -- No functions are available yet.

    -- * Auxiliary functions
    product',
  )
where

import Numeric.Log
import qualified Statistics.Distribution as S
import qualified Statistics.Distribution.Exponential as S
import qualified Statistics.Distribution.Gamma as S
import qualified Statistics.Distribution.Normal as S

-- | Improper uniform probability density function; larger than 0.
positive :: Double -> Log Double
positive x | x <= 0 = 0
           | otherwise = 1

-- | Improper uniform probability density function; lower than 0.
negative :: Double -> Log Double
negative x | x >= 0 = 0
           | otherwise = 1

-- | Uniform probability density in [a, b].
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

-- | Normal distributed density.
normal ::
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  Double ->
  Log Double
normal m s x = Exp $ S.logDensity d x
  where d = S.normalDistr m s

-- | Exponential distributed density.
exponential ::
  -- | Rate.
  Double ->
  Double ->
  Log Double
exponential l x = Exp $ S.logDensity d x
  where
    d = S.exponential l

-- | Gamma distributed density.
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

-- | Intelligent product that stops when encountering a zero.
product' :: [Log Double] -> Log Double
product' [] = 0
product' [x] = x
product' (x : xs)
  | x == 0 = 0
  | otherwise = x * product' xs
