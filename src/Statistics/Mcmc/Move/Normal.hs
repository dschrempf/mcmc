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
  , moveNormalDouble
  ) where

import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC

import Statistics.Mcmc.Types

{-# INLINE delta #-}
delta :: (a -> Double) -> (Double -> a -> a) -> NormalDistribution -> a -> GenIO -> IO a
delta get set d x g = do
  dx <- genContinuous d g
  return $ set (get x + dx) x

{-# INLINE logDens #-}
logDens :: (a -> Double) -> NormalDistribution -> a -> a -> Log Double
logDens get d x y = Exp $ logDensity d (get y - get x)

-- | A symmetric move with normally distributed density.
moveNormal
  :: (a -> Double)      -- ^ Getter.
  -> (Double -> a -> a) -- ^ Setter.
  -> String             -- ^ Name.
  -> Double             -- ^ Mean.
  -> Double             -- ^ Standard deviation.
  -> Move a
moveNormal get set n m s = Move n (delta get set d) (logDens get d)
  where d = normalDistr m s

-- | A symmetric move with normally distributed density; specialized to a one
-- dimensional state space of type 'Double'; see 'moveNormal'.
moveNormalDouble :: String -> Double -> Double -> Move Double
moveNormalDouble = moveNormal id const
