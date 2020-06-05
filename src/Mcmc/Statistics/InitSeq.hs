{-# LANGUAGE FlexibleContexts #-}

{- |
Module      :  Mcmc.Statistics.InitSeq
Description :  Initial convex sequence estimator of asymptotic variance
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed Jun  3 18:38:18 2020.

See Geyer, C. J., Introduction to Markov chain Monte Carlo, In Handbook of
Markov Chain Monte Carlo (pp. 45) (2011). Chapman \& Hall/CRC.

-}

-- XXX: All this is slow:
--
-- - snoc is used extensively
--
-- - vectors are manipulated extensively
--
-- - functions recompute the length repeatedly
--
-- Solution: Use mutable vectors and store indices.

module Mcmc.Statistics.InitSeq
  ( initSeq
  ) where

import qualified Data.Vector.Generic as V
import Data.Vector.Generic (Vector)
import Statistics.Autocorrelation

-- Sum of two consecutive elements.
sum2 :: Vector v Double => Int -> v Double -> Double
sum2 i = V.sum . V.slice i 2

-- Takes the vector of autocovariances and computes the vector Gamma_k; Eq.
-- (1.18), page 16.
gamma :: (Vector v Double, Vector v Int) => v Double -> v Double
gamma c = V.generate l (`sum2` c)
  where l = (V.length c + 1) `quot` 2

-- Truncated Gamma function; not numbered.
gammaTruncated :: (Vector v Double, Vector v Int) => v Double -> v Double
gammaTruncated c = V.takeWhile (>0) (gamma c) `V.snoc` 0

-- Differences between values in vector.
diff :: (Num a, Vector v a) => v a -> v a
diff v = V.zipWith (-) (V.tail v) v
{-# SPECIALIZE diff :: (Vector v Double) => v Double -> v Double #-}

-- Slope of last element.
slopeLast :: (Fractional a, Vector v a) => v a -> v a -> a
slopeLast xs ys = (V.last ys - secondLast ys) / (V.last xs - secondLast xs)
{-# SPECIALIZE slopeLast :: (Vector v Double) => v Double -> v Double -> Double #-}

secondLast :: Vector v a => v a -> a
secondLast = V.last . V.init
{-# SPECIALIZE secondLast :: (Vector v Double) => v Double -> Double #-}

rmSecondLast :: Vector v a => v a -> v a
rmSecondLast v = V.slice 0 (l-2) v `V.snoc` V.last v
  where l = V.length v
{-# SPECIALIZE rmSecondLast :: (Vector v Double) => v Double -> v Double #-}

-- Pool the last value in a vector until convexity is preserved.
pool :: (Fractional a, Ord a, Vector v a) => v a -> v a -> (v a, v a)
pool xs ys = if s0 > s1 then (xs, ys) else pool (rmSecondLast xs) (rmSecondLast ys)
  where
    s0 = slopeLast xs ys
    s1 = slopeLast (V.init xs) (V.init ys)
{-# SPECIALIZE pool :: (Vector v Double) => v Double -> v Double -> (v Double, v Double) #-}

-- Greatest convex minorant of a list of observations. Uses the Pool Adjacent
-- Violators Algorithm (PAVA). The predictors are set to
--
-- \[
-- x_0 = 0, x_1 = 1, \ldots, x_{l-1} = l-1,
-- \]
--
-- where $l$ is the length of the given observation vector. Hence, the
-- predictors are ordered and there are no ties, i.e.,
--
-- \[
-- x_0 < x_1 < ... < x_{l-1}.
-- \]
--
-- See https://en.wikipedia.org/wiki/Isotonic_regression, although this function
-- is antitonic in that the resulting vector is decreasing.
--
-- XXX: Especially here, mutable vectors should be used.
--
-- TODO: Consistent naming of variables:
--
-- - xs and ys: x and y values of gcm
--
-- - zs and os: x and y values of observations (predictors and responses)
gcm :: (Vector v Double, Vector v Int) => v Double -> v Double
gcm os = go (V.singleton 0 :: v Int) (V.head os) 0 (V.tail os)
  where
    l  = V.length os
    go xs ys i os = undefined

-- | Initial convex sequence estimator of the asymptotic variance. Eq. (1.19),
-- page 16.
initSeq :: (Vector v Double, Vector v Int) => v Double -> Double
initSeq v = V.head c + 2 * V.sum (gcm . gammaTruncated $ c)
  where c = autocovariance v
