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

import qualified Data.Vector.Generic as V
import Data.Vector.Generic (Vector)

import Statistics.Gcm
import Statistics.Autocorrelation

-- TODO: Test and check these functions.

-- Sum of two consecutive elements.
sum2 :: Vector v Double => Int -> v Double -> Double
sum2 i = V.sum . V.slice (2*i) 2

-- Takes the vector of autocovariances and computes the vector Gamma_k; Eq.
-- (1.18), page 16.
gamma :: (Vector v Double, Vector v Int) => v Double -> v Double
gamma c = V.generate l (`sum2` c)
  -- TODO: +1?
  where l = (V.length c + 1) `quot` 2

-- TODO: Use mutable vector so that not the complete gamma vector needs to be created.

-- Truncated Gamma function; not numbered.
gammaTruncated :: (Vector v Double, Vector v Int) => v Double -> v Double
gammaTruncated c = V.takeWhile (>0) (gamma c) `V.snoc` 0

-- | Initial convex sequence estimator of the asymptotic variance. Eq. (1.19),
-- page 16.
initSeq :: (Vector v Double, Vector v Int) => v Double -> Double
initSeq v = V.head c + 2 * V.sum (uncurry smooth . gcm . gammaTruncated $ c)
  where c = autocovariance v
