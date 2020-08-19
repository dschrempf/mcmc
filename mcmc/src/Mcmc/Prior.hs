{-# LANGUAGE BangPatterns #-}

-- |
-- Module      :  Prior
-- Description :  Convenience functions to compute priors
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 13:26:14 2020.
module Mcmc.Prior
  ( -- * Continuous priors
    positive,
    negative,
    uniform,
    normal,
    exponential,
    gamma,

    -- * Discrete priors

    -- No discrete priors are available yet.

    -- * Auxiliary functions
    product',
  )
where

import Control.Monad
import Data.Maybe (fromMaybe)
import Numeric.Log
import qualified Statistics.Distribution as S
import qualified Statistics.Distribution.Exponential as S
import qualified Statistics.Distribution.Gamma as S
import qualified Statistics.Distribution.Normal as S

-- | Improper uniform prior; larger than 0.
positive :: Double -> Log Double
positive x
  | x <= 0 = 0
  | otherwise = 1

-- | Improper uniform prior; lower than 0.
negative :: Double -> Log Double
negative x
  | x >= 0 = 0
  | otherwise = 1

-- | Uniform prior on [a, b].
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
  where
    d = S.normalDistr m s

-- | Exponential distributed prior.
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

-- | Intelligent product that stops when encountering a zero.
--
-- Use with care because the elements have to be checked for positiveness, and
-- this can take some time if the list is long and does not contain any zeroes.
product' :: [Log Double] -> Log Double
product' = fromMaybe 0 . prodM

-- The type could be generalized to any MonadPlus Integer
prodM :: [Log Double] -> Maybe (Log Double)
prodM = foldM (\ !acc x -> (acc * x) <$ guard (acc /= 0)) 1
