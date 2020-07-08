{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal.Bactrian
-- Description :  Bactrian proposals
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
module Mcmc.Proposal.Bactrian
  ( slideBactrian,
    scaleBactrian,
  )
where

import Lens.Micro
import Mcmc.Proposal
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
  ProposalSimple a
bactrianAdditiveSimple l m s t
  | m < 0 = error "bactrianAdditiveSimple: Spike parameter negative."
  | m >= 1 = error "bactrianAdditiveSimple: Spike parameter 1.0 or larger."
  | s <= 0 = error "bactrianAdditiveSimple: Standard deviation 0.0 or smaller."
  | otherwise = ProposalSimple (bactrianAdditive l m (t * s)) Nothing

-- | Additive symmetric proposal with kernel similar to the silhouette of a
-- Bactrian camel. The [Bactrian
-- kernel](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/figure/fig01)
-- is a mixture of two symmetrically arranged normal distributions. The spike
-- parameter loosely determines the standard deviations of the individual humps
-- while the other parameter refers to the standard deviation of the complete
-- Bactrian kernel.
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
  Proposal a
slideBactrian n w l m s t = Proposal n w (bactrianAdditiveSimple l m s 1.0) tnr
  where
    tnr = if t then Just $ tuner (bactrianAdditiveSimple l m s) else Nothing

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

bactrianKernel :: Double -> Double -> Double -> Log Double
bactrianKernel m s x = Exp $ log $ kernel1 + kernel2
  where
    mn = m * s
    sd = sqrt (1 - m * m) * s
    dist1 = normalDistr (- mn) sd
    dist2 = normalDistr mn sd
    kernel1 = density dist1 x
    kernel2 = density dist2 x

bactrianMultKernel :: Lens' a Double -> Double -> Double -> a -> a -> Log Double
bactrianMultKernel l m s x y = bactrianKernel m s dx
  where
    dx = y ^. l / x ^. l

bactrianMultSimple :: Lens' a Double -> Double -> Double -> Double -> ProposalSimple a
bactrianMultSimple l m s t
  | m < 0 = error "bactrianMultSimple: Spike parameter negative."
  | m >= 1 = error "bactrianMultSimple: Spike parameter 1.0 or larger."
  | s <= 0 = error "bactrianMultSimple: Standard deviation 0.0 or smaller."
  | otherwise = ProposalSimple (bactrianMult l m (t * s)) (Just $ bactrianMultKernel l m (t * s))

-- | Multiplicative proposal with kernel similar to the silhouette of a Bactrian
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
  Proposal a
scaleBactrian n w l m s t = Proposal n w (bactrianMultSimple l m s 1.0) tnr
  where
    tnr = if t then Just $ tuner (bactrianMultSimple l m s) else Nothing
