{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

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

module Mcmc.Statistics.InitSeq
  (
    initSeq
  ) where

import qualified Data.Vector.Unboxed as V
import Data.Vector.Unboxed (Vector, Unbox)

import Statistics.Autocorrelation
import Statistics.Gcm
import Statistics.Pava.Common

-- TODO: Test and check these functions.

-- Sum of two consecutive elements.
sum2 :: Int -> Vector Double -> Double
sum2 i = V.sum . V.slice (2*i) 2

-- Takes the vector of autocovariances and computes the vector Gamma_k; Eq.
-- (1.18), page 16.
gamma :: Vector Double -> Vector Double
gamma c = V.generate l (`sum2` c)
  where l = V.length c `quot` 2

-- Truncated Gamma function; not numbered.
gammaTruncated :: Vector Double -> Vector Double
-- XXX: Use mutable vector so that not the complete gamma vector needs to be created.
gammaTruncated c = V.takeWhile (>0) (gamma c) `V.snoc` 0

-- | Initial convex sequence estimator of the asymptotic variance. Eq. (1.19),
-- page 16.
initSeq :: (Real a, Unbox a) => Vector a -> Double
initSeq v = (- V.head c) + 2 * V.sum (smooth (V.fromList gcmIs) (V.fromList gcmVals))
  where c = autocovariance (V.map realToFrac v)
        gs = gammaTruncated c
        is = V.fromList [1 .. V.length gs]
        (gcmIs, gcmVals, _) = gcm is gs
