{-# LANGUAGE RankNTypes #-}

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

import Lens.Micro
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC

import Statistics.Mcmc.Types

delta :: Lens' a Double -> NormalDistribution -> a -> GenIO -> IO a
delta l d x g = do
  dx <- genContVar d g
  return $ set l (x^.l + dx) x
{-# INLINEABLE delta #-}

-- XXX: Technically, only a Getter is needed here.
logDens :: Lens' a Double -> NormalDistribution -> a -> a -> Log Double
logDens l d x y = Exp $ logDensity d (y^.l - x^.l)
{-# INLINEABLE logDens #-}

-- | A symmetric move with normally distributed density.
moveNormal
  :: Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> Double             -- ^ Mean.
  -> Double             -- ^ Standard deviation.
  -> Move a
moveNormal l n m s = Move n (delta l d) (logDens l d)
  where d = normalDistr m s

-- Ridiculous.
lid :: Lens' a a
lid f = f

-- | A symmetric move with normally distributed density; specialized to a one
-- dimensional state space of type 'Double'; see 'moveNormal'.
moveNormalDouble :: String -> Double -> Double -> Move Double
moveNormalDouble = moveNormal lid

-- -- It turns out that direct implementation is not faster (or only marginally
-- -- faster); tested with criterion for 50000 jumps.

-- deltaDouble :: NormalDistribution -> Double -> Double -> GenIO -> IO Double
-- deltaDouble d m x g = do
--   dx <- genContVar d g
--   return $ (x + m + dx)
-- {-# INLINEABLE deltaDouble #-}

-- logDensDouble :: NormalDistribution -> Double -> Double -> Log Double
-- logDensDouble d x y = Exp $ logDensity d (y - x)
-- {-# INLINEABLE logDensDouble #-}

-- -- | A symmetric move with normally distributed density.
-- moveNormalDoubleFast
--   :: String             -- ^ Name.
--   -> Double             -- ^ Mean.
--   -> Double             -- ^ Standard deviation.
--   -> Move Double
-- moveNormalDoubleFast n m s = Move n (deltaDouble d m) (logDensDouble d)
--   where d = normalDistr 0 s
