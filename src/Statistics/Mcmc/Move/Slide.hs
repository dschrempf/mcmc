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
  , slideUniform
  , slideUniformDouble
  ) where

import Lens.Micro
import Statistics.Distribution.Normal
import Statistics.Distribution.Uniform

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
  -> Int                -- ^ Weight.
  -> Lens' a Double     -- ^ Instruction about which parameter to change.
  -> Double             -- ^ Mean.
  -> Double             -- ^ Standard deviation.
  -> Bool               -- ^ Enable tuning.
  -> Move a
slide n w l m s True  = Move n w (slideSimple l m s 1.0) (Just $ tuner (slideSimple l m s))
slide n w l m s False = Move n w (slideSimple l m s 1.0)  Nothing

-- | See 'slide'; specialized to a one dimensional state space of type 'Double'.
slideDouble
  :: String -- ^ Name.
  -> Int    -- ^ Weight.
  -> Double -- ^ Mean.
  -> Double -- ^ Standard deviation.
  -> Bool   -- ^ Enable tuning.
  -> Move Double
slideDouble n w = slide n w id

-- The actual move with tuning parameter.
slideUniformSimple
  :: Lens' a Double
  -> Double
  -> Double
  -> MoveSimple a
slideUniformSimple l d t = moveGenericContinuous l (uniformDistr (-t*d) (t*d)) (+) (-)

-- | Additive move with uniformly distributed density.
slideUniform
  :: String             -- ^ Name.
  -> Int                -- ^ Weight.
  -> Lens' a Double     -- ^ Instruction about which parameter to change.
  -> Double             -- ^ Delta.
  -> Bool               -- ^ Enable tuning.
  -> Move a
slideUniform n w l d True  = Move n w (slideUniformSimple l d 1.0)
                           (Just $ tuner (slideUniformSimple l d))
slideUniform n w l d False = Move n w (slideUniformSimple l d 1.0)  Nothing

-- | See 'slideUniform'; specialized to a one dimensional state space of type
-- 'Double'.
slideUniformDouble
  :: String -- ^ Name.
  -> Int    -- ^ Weight.
  -> Double -- ^ Delta.
  -> Bool   -- ^ Enable tuning.
  -> Move Double
slideUniformDouble n w = slideUniform n w id
