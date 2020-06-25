{-# LANGUAGE RankNTypes #-}

-- TODO: Test this move!

-- TODO:
-- slideBactrian
--
-- scaleBactrian

-- |
-- Module      :  Mcmc.Move.Bactrian
-- Description :  Bactrian moves
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jun 25 15:49:48 2020.
--
-- See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/.
module Mcmc.Move.Bactrian
  ( bactrian
  )
where

import Lens.Micro
import Mcmc.Move
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC
import System.Random.MWC.Distributions

-- bactrianSample lens spike stdDev
bactrianSample :: Lens' a Double -> Double -> Double -> a -> GenIO -> IO a
bactrianSample l m s x g = do
  let v = x ^. l
      d = normalDistr (m * s) (sqrt (1 - m * m) * s)
  dv <- genContVar d g
  b <- bernoulli 0.5 g
  let v' = if b then v + dv else v - dv
  return $ x & l .~ v'

-- bactrianSimple lens spike stdDev tune
bactrianSimple :: Lens' a Double -> Double -> Double -> Double -> MoveSimple a
bactrianSimple l m s t | m < 0 = error "bactrianSimple: Spike parameter negative."
                       | m >= 1 = error "bactrianSimple: Spike parameter 1.0 or larger."
                       | s <= 0 = error "bactrianSimple: Standard deviation 0.0 or smaller."
                       | otherwise = MoveSimple (bactrianSample l m (t*s)) Nothing

-- TODO: Improve help.
-- | Additive move with density similar to a Bactrian camel. See
-- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/.
bactrian ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Spike parameter.
  Double ->
  -- | Standard deviation.
  Double ->
  -- TODO: Check how RevBayes tunes the Bactrian move.
  -- | Enable tuning.
  Bool ->
  Move a
-- TODO: This and all other move definitions can be condensed. Only the tuner
-- has to be changed according to the given Bool parameter.
bactrian n w l m s True =
  Move n w (bactrianSimple l m s 1.0) (Just $ tuner (bactrianSimple l m s))
bactrian n w l m s False = Move n w (bactrianSimple l m s 1.0) Nothing
