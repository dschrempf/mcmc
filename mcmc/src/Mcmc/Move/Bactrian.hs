{-# LANGUAGE RankNTypes #-}

-- TODO: Test these moves!

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
    scaleBactrian,
  )
where

import Lens.Micro
import Mcmc.Move
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC
import System.Random.MWC.Distributions

bactrianSample ::
  Double ->
  Double ->
  GenIO ->
  IO Double
bactrianSample m s g = do
  let mn = m * s
      sd = sqrt (1 - m * m) * s
      d = normalDistr mn sd
  x <- genContVar d g
  b <- bernoulli 0.5 g
  return $ if b then x else (- x)

bactrianAdditive ::
  Lens' a Double ->
  Double ->
  Double ->
  a ->
  GenIO ->
  IO a
bactrianAdditive l m s x g = do
  dx <- bactrianSample m s g
  return $ x & l +~ dx

-- bactrianSimple lens spike stdDev tune forwardOp backwardOp
bactrianAdditiveSimple ::
  Lens' a Double ->
  Double ->
  Double ->
  Double ->
  MoveSimple a
bactrianAdditiveSimple l m s t
  | m < 0 = error "bactrianAdditiveSimple: Spike parameter negative."
  | m >= 1 = error "bactrianAdditiveSimple: Spike parameter 1.0 or larger."
  | s <= 0 = error "bactrianAdditiveSimple: Standard deviation 0.0 or smaller."
  | otherwise = MoveSimple (bactrianAdditive l m (t * s)) Nothing

-- | Additive symmetric move with density similar to the silhouette of a
-- Bactrian camel. The [Bactrian
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
  Move n w (bactrianAdditiveSimple l m s 1.0) (Just $ tuner (bactrianAdditiveSimple l m s))
slideBactrian n w l m s False = Move n w (bactrianAdditiveSimple l m s 1.0) Nothing

bactrianMult ::
  Lens' a Double ->
  Double ->
  Double ->
  a ->
  GenIO ->
  IO a
bactrianMult l m s x g = do
  dx <- bactrianSample m s g
  return $ x & l %~ (* (1 + dx))

bactrianDens :: Double -> Double -> Double -> Log Double
bactrianDens m s x = Exp $ log $ dens1 + dens2
  where
    mn = m * s
    sd = sqrt (1 - m * m) * s
    dist1 = normalDistr (- mn) sd
    dist2 = normalDistr mn sd
    dens1 = density dist1 x
    dens2 = density dist2 x

bactrianMultDens :: Lens' a Double -> Double -> Double -> a -> a -> Log Double
bactrianMultDens l m s x y = bactrianDens m s dx
  where
    dx = y ^. l / x ^. l

bactrianMultSimple :: Lens' a Double -> Double -> Double -> Double -> MoveSimple a
bactrianMultSimple l m s t
  | m < 0 = error "bactrianMultSimple: Spike parameter negative."
  | m >= 1 = error "bactrianMultSimple: Spike parameter 1.0 or larger."
  | s <= 0 = error "bactrianMultSimple: Standard deviation 0.0 or smaller."
  | otherwise = MoveSimple (bactrianMult l m (t * s)) (Just $ bactrianMultDens l m (t * s))

-- | Multiplicative move with density similar to the silhouette of a Bactrian
-- camel. See 'slideBactrian'.
scaleBactrian ::
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
scaleBactrian n w l m s True =
  Move n w (bactrianMultSimple l m s 1.0) (Just $ tuner (bactrianMultSimple l m s))
scaleBactrian n w l m s False = Move n w (bactrianMultSimple l m s 1.0) Nothing
