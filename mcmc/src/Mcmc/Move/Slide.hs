{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Move.Slide
-- Description :  Normally distributed move
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May  6 10:59:13 2020.
module Mcmc.Move.Slide
  ( slide,
    slideDouble,
    slideStandard,
    slideUniform,
    slideUniformDouble,
  )
where

import Lens.Micro
import Mcmc.Move
import Mcmc.Move.Generic
import Statistics.Distribution.Normal
import Statistics.Distribution.Uniform

-- The actual move with tuning parameter.
slideSimple :: Lens' a Double -> Double -> Double -> Double -> MoveSimple a
slideSimple l m s t = moveGenericContinuous l (normalDistr m (s * t)) (+) (-)

-- | Additive move with normally distributed density.
slide ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Move a
slide n w l m s True =
  Move n w (slideSimple l m s 1.0) (Just $ tuner (slideSimple l m s))
slide n w l m s False = Move n w (slideSimple l m s 1.0) Nothing

-- | See 'slide'; specialized to a one dimensional state space of type 'Double'.
slideDouble ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Move Double
slideDouble n w = slide n w id

-- | Additive move with normally distributed density with a mean of 0.0 and a
-- standard deviation of 1.0.
slideStandard ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Enable tuning.
  Bool ->
  Move a
slideStandard n w l True =
  Move n w (slideSimple l 0.0 1.0 1.0) (Just $ tuner (slideSimple l 0.0 1.0))
slideStandard n w l False = Move n w (slideSimple l 0.0 1.0 1.0) Nothing

-- The actual move with tuning parameter.
slideUniformSimple :: Lens' a Double -> Double -> Double -> MoveSimple a
slideUniformSimple l d t =
  moveGenericContinuous l (uniformDistr (- t * d) (t * d)) (+) (-)

-- | Additive move with uniformly distributed density.
slideUniform ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Delta.
  Double ->
  -- | Enable tuning.
  Bool ->
  Move a
slideUniform n w l d True =
  Move n w (slideUniformSimple l d 1.0) (Just $ tuner (slideUniformSimple l d))
slideUniform n w l d False = Move n w (slideUniformSimple l d 1.0) Nothing

-- | See 'slideUniform'; specialized to a one dimensional state space of type
-- 'Double'.
slideUniformDouble ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Delta.
  Double ->
  -- | Enable tuning.
  Bool ->
  Move Double
slideUniformDouble n w = slideUniform n w id
