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
    slideSymmetric,
    slideUniform,
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
slide n w l m s t = Move n w (slideSimple l m s 1.0) tnr
  where tnr = if t then  Just $ tuner (slideSimple l m s) else Nothing

-- The actual move with tuning parameter.
slideSymmetricSimple :: Lens' a Double -> Double -> Double -> MoveSimple a
slideSymmetricSimple l s t = moveSymmetricGenericContinuous l (normalDistr 0.0 (s * t)) (+)

-- | Additive move with normally distributed density with mean zero. This move
-- is very fast, because the Metropolis-Hastings ratio does not include
-- calculation of the forwards and backwards densities.
slideSymmetric ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Move a
slideSymmetric n w l s t = Move n w (slideSymmetricSimple l s 1.0) tnr
  where tnr = if t then Just $ tuner (slideSymmetricSimple l s) else Nothing

-- The actual move with tuning parameter.
slideUniformSimple :: Lens' a Double -> Double -> Double -> MoveSimple a
slideUniformSimple l d t =
  moveSymmetricGenericContinuous l (uniformDistr (- t * d) (t * d)) (+)

-- | Additive move with uniformly distributed density. This move is very fast,
-- because the Metropolis-Hastings ratio does not include calculation of the
-- forwards and backwards densities.
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
slideUniform n w l d t = Move n w (slideUniformSimple l d 1.0) tnr
  where tnr = if t then Just $ tuner (slideUniformSimple l d) else Nothing
