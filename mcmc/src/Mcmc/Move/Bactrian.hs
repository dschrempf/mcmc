{-# LANGUAGE RankNTypes #-}

-- TODO: Test this move!

-- TODO: ScaleBactrian.

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
  ( slideBactrian,
  )
where

import Lens.Micro
import Mcmc.Move
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC
import System.Random.MWC.Distributions

-- bactrianSample lens spike stdDev forwardOp backwardOp
bactrianSample ::
  Lens' a Double ->
  Double ->
  Double ->
  (Double -> Double -> Double) ->
  (Double -> Double -> Double) ->
  a ->
  GenIO ->
  IO a
bactrianSample l m s fop bop x g = do
  let v = x ^. l
      d = normalDistr (m * s) (sqrt (1 - m * m) * s)
  dv <- genContVar d g
  b <- bernoulli 0.5 g
  let v' = if b then fop v dv else bop v dv
  return $ x & l .~ v'

-- bactrianSimple lens spike stdDev tune forwardOp backwardOp
bactrianSimple ::
  Lens' a Double ->
  Double ->
  Double ->
  (Double -> Double -> Double) ->
  (Double -> Double -> Double) ->
  Double ->
  MoveSimple a
bactrianSimple l m s fop bop t
  | m < 0 = error "bactrianSimple: Spike parameter negative."
  | m >= 1 = error "bactrianSimple: Spike parameter 1.0 or larger."
  | s <= 0 = error "bactrianSimple: Standard deviation 0.0 or smaller."
  | otherwise = MoveSimple (bactrianSample l m (t * s) fop bop) Nothing

-- | Additive move with density similar to the silhouette of a Bactrian camel.
-- The [Bactrian
-- density](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/figure/fig01)
-- is a mixture of two symmetrically arranged normal distributions. The spike
-- parameter loosely determines the standard deviations of the individual humps
-- while the other parameter refers to the standard deviation of the complete
-- Bactrian density.
--
-- See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/.
slideBactrian ::
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
  -- | Enable tuning.
  Bool ->
  Move a
-- TODO: This and all other move definitions can be condensed. Only the tuner
-- has to be changed according to the given Bool parameter.
slideBactrian n w l m s True =
  Move n w (bactrianSimple l m s (+) (-) 1.0) (Just $ tuner (bactrianSimple l m s (+) (-)))
slideBactrian n w l m s False = Move n w (bactrianSimple l m s (+) (-) 1.0) Nothing
