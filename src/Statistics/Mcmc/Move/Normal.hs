{- |
Module      :  Statistics.Mcmc.Move.Normal
Description :  Normally distributed move
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 10:59:13 2020.

-}

module Statistics.Mcmc.Move.Normal
  ( moveNormal
  ) where

import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC

import Statistics.Mcmc.Types

{-# INLINE delta #-}
delta :: NormalDistribution -> Double -> GenIO -> IO Double
delta d x g = do
  dx <- genContinuous d g
  return $ x + dx

{-# INLINE logDens #-}
logDens :: NormalDistribution -> Double -> Double -> Log Double
logDens d x y = Exp $ logDensity d (y - x)

-- | A symmetric move with normally distributed density.
moveNormal
  :: Double -- ^ Mean.
  -> Double -- ^ Standard deviation.
  -> Move Double
moveNormal m s = Move (delta d) (logDens d)
  where d = normalDistr m s
