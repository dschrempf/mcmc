{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Statistics.Mcmc.Moves.Generic
Description :  Generic interface to create moves
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 14 20:26:27 2020.

-}

module Statistics.Mcmc.Moves.Generic
  ( moveGeneric
  ) where

import Lens.Micro
import Numeric.Log
import Statistics.Distribution
import System.Random.MWC

import Statistics.Mcmc.Types

jump
  :: (ContDistr d, ContGen d)
  => Lens' a Double
  -> d
  -> (Double -> Double -> Double)
  -> a
  -> GenIO
  -> IO a
jump l d f x g = do
  dx <- genContVar d g
  return $ set l ((x^.l) `f` dx) x
{-# INLINEABLE jump #-}

-- XXX: Technically, only a Getter is needed here.
logDens
  :: (ContDistr d, ContGen d)
  => Lens' a Double
  -> d
  -> (Double -> Double -> Double)
  -> a -> a -> Log Double
logDens l d fInv x y = Exp $ logDensity d ((y^.l) `fInv` (x^.l))
{-# INLINEABLE logDens #-}

-- | A symmetric move with normally distributed density.
moveGeneric
  :: (ContDistr d, ContGen d)
  => Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> d                  -- ^ Probability distribution.
  -> (Double -> Double -> Double) -- ^ Forward operator, e.g. (+), so that x + dx = y.
  -> (Double -> Double -> Double) -- ^ Inverse operator, e.g.,(-), so that y - dx = x.
  -> Move a
moveGeneric l n d f fInv = Move n (jump l d f) (logDens l d fInv)
