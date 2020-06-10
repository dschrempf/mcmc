{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Mcmc.Move.Generic
Description :  Generic interface to create moves
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 14 20:26:27 2020.

-}

module Mcmc.Move.Generic
  ( moveGenericContinuous
  -- , moveGenericDiscrete
  )
where

import           Numeric.Log
import           Statistics.Distribution
import           System.Random.MWC

import           Mcmc.Move

jumpCont
  :: (ContDistr d, ContGen d)
  => (Double -> a -> a)
  -> (a -> Double)
  -> d
  -> (Double -> Double -> Double)
  -> a
  -> GenIO
  -> IO a
jumpCont setter getter d f x g = do
  dx <- genContVar d g
  return $ setter (getter x `f` dx) x
{-# INLINABLE jumpCont #-}

logDensCont
  :: (ContDistr d, ContGen d)
  => (a -> Double)
  -> d
  -> (Double -> Double -> Double)
  -> a
  -> a
  -> Log Double
logDensCont getter d fInv x y = Exp $ logDensity d (getter y `fInv` getter x)
{-# INLINABLE logDensCont #-}

-- | Generic function to create moves for continuous parameters ('Double').
moveGenericContinuous
  :: (ContDistr d, ContGen d)
  => (Double -> a -> a)           -- ^ Setter.
  -> (a -> Double)                -- ^ Getter.
  -> d                            -- ^ Probability distribution
  -> (Double -> Double -> Double) -- ^ Forward operator, e.g. (+), so that x + dx = y.
  -> (Double -> Double -> Double) -- ^ Inverse operator, e.g.,(-), so that y - dx = x.
  -> MoveSimple a
moveGenericContinuous setter getter d f fInv =
  MoveSimple (jumpCont setter getter d f) (logDensCont getter d fInv)

-- jumpDiscrete
--   :: (DiscreteDistr d, DiscreteGen d)
--   => Lens' a Int
--   -> d
--   -> (Int -> Int -> Int)
--   -> a
--   -> GenIO
--   -> IO a
-- jumpDiscrete l d f x g = do
--   dx <- genDiscreteVar d g
--   return $ set l ((x ^. l) `f` dx) x
-- {-# INLINABLE jumpDiscrete #-}

-- -- XXX: Technically, only a Getter is needed here.
-- logDensDiscrete
--   :: (DiscreteDistr d, DiscreteGen d)
--   => Lens' a Int
--   -> d
--   -> (Int -> Int -> Int)
--   -> a
--   -> a
--   -> Log Double
-- logDensDiscrete l d fInv x y =
--   Exp $ logProbability d ((y ^. l) `fInv` (x ^. l))
-- {-# INLINABLE logDensDiscrete #-}

-- -- | Generic function to create moves for discrete parameters ('Int').
-- moveGenericDiscrete
--   :: (DiscreteDistr d, DiscreteGen d)
--   => Lens' a Int         -- ^ Instruction about which parameter to change.
--   -> d                   -- ^ Probability distribution.
--   -> (Int -> Int -> Int) -- ^ Forward operator, e.g. (+), so that x + dx = y.
--   -> (Int -> Int -> Int) -- ^ Inverse operator, e.g.,(-), so that y - dx = x.
--   -> MoveSimple a
-- moveGenericDiscrete l fd f fInv =
--   MoveSimple (jumpDiscrete l fd f) (logDensDiscrete l fd fInv)
