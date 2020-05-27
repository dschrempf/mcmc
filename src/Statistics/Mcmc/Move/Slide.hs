{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Statistics.Mcmc.Move.Slide
Description :  Normally distributed move
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 10:59:13 2020.

-}

module Statistics.Mcmc.Move.Slide
  ( slide
  , slideDouble
  ) where

import Lens.Micro
import Statistics.Distribution.Normal

import Statistics.Mcmc.Move.Generic
import Statistics.Mcmc.Move

-- The actual move with tuning parameter.
slideSimple
  :: Lens' a Double
  -> Double
  -> Double
  -> Double
  -> MoveSimple a
slideSimple l m s t = moveGenericContinuous l (normalDistr m (s*t)) (+) (-)

-- | Additive move with normally distributed density.
slide
  :: String             -- ^ Name.
  -> Lens' a Double     -- ^ Instruction about which parameter to change.
  -> Double             -- ^ Mean.
  -> Double             -- ^ Standard deviation.
  -> Bool               -- ^ Enable tuning.
  -> Move a
slide n l m s True  = Move n (slideSimple l m s 1.0) (Just $ tuner (slideSimple l m s))
slide n l m s False = Move n (slideSimple l m s 1.0)  Nothing

-- | Additive move with normally distributed density; specialized to a one
-- dimensional state space of type 'Double'.
slideDouble
  :: String -- ^ Name.
  -> Double -- ^ Mean.
  -> Double -- ^ Standard deviation.
  -> Bool   -- ^ Enable tuning.
  -> Move Double
slideDouble n = slide n id
