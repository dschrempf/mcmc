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
  ( SpikeParameter,
    slideBactrian,
    scaleBactrian,
  )
where

import Mcmc.Proposal
import Mcmc.Statistics.Types
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC
import System.Random.MWC.Distributions

-- | Type synonym indicating the spike parameter.
type SpikeParameter = Double

genBactrian ::
  SpikeParameter ->
  StandardDeviation ->
  GenIO ->
  IO Double
genBactrian m s g = do
  let mn = m * s
      sd = sqrt (1 - m * m) * s
      d = normalDistr mn sd
  x <- genContVar d g
  b <- bernoulli 0.5 g
  return $ if b then x else - x

logDensityBactrian :: SpikeParameter -> StandardDeviation  -> Double -> Log Double
logDensityBactrian m s x = Exp $ log $ kernel1 + kernel2
  where
    mn = m * s
    sd = sqrt (1 - m * m) * s
    dist1 = normalDistr (- mn) sd
    dist2 = normalDistr mn sd
    kernel1 = density dist1 x
    kernel2 = density dist2 x

bactrianAdditive ::
  SpikeParameter ->
  StandardDeviation  ->
  ProposalSimple Double
bactrianAdditive m s x g = do
  dx <- genBactrian m s g
  return (x + dx, 1.0, 1.0)

-- bactrianSimple lens spike stdDev tune forwardOp backwardOp
bactrianAdditiveSimple ::
  SpikeParameter ->
  StandardDeviation  ->
  TuningParameter  ->
  ProposalSimple Double
bactrianAdditiveSimple m s t
  | m < 0 = error "bactrianAdditiveSimple: Spike parameter negative."
  | m >= 1 = error "bactrianAdditiveSimple: Spike parameter 1.0 or larger."
  | s <= 0 = error "bactrianAdditiveSimple: Standard deviation 0.0 or smaller."
  | otherwise = bactrianAdditive m (t * s)

-- | Additive symmetric proposal with kernel similar to the silhouette of a
-- Bactrian camel.
--
-- The [Bactrian
-- kernel](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/figure/fig01) is
-- a mixture of two symmetrically arranged normal distributions. The spike
-- parameter \(m \in (0, 1)\) loosely determines the standard deviations of the
-- individual humps while the second parameter \(s > 0\) refers to the
-- standard deviation of the complete Bactrian kernel.
--
-- See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/.
slideBactrian ::
  SpikeParameter ->
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal Double
slideBactrian m s = createProposal description (bactrianAdditiveSimple m s) (PDimension 1)
  where
    description = PDescription $ "Slide Bactrian; spike: " ++ show m ++ ", sd: " ++ show s

-- We have:
-- x  (1+dx ) = x'
-- x' (1+dx') = x.
--
-- Hence,
-- dx' = 1/(1-dx) - 1.
fInv :: Double -> Double
fInv dx = recip (1 - dx) - 1

bactrianMult ::
  SpikeParameter ->
  StandardDeviation  ->
  ProposalSimple Double
bactrianMult m s x g = do
  du <- genBactrian m s g
  let qXY = logDensityBactrian m s du
      qYX = logDensityBactrian m s (fInv du)
      u = 1.0 + du
      jac = Exp $ log $ recip u
  return (x * u, qYX / qXY, jac)

bactrianMultSimple :: SpikeParameter -> StandardDeviation  -> TuningParameter  -> ProposalSimple Double
bactrianMultSimple m s t
  | m < 0 = error "bactrianMultSimple: Spike parameter negative."
  | m >= 1 = error "bactrianMultSimple: Spike parameter 1.0 or larger."
  | s <= 0 = error "bactrianMultSimple: Standard deviation 0.0 or smaller."
  | otherwise = bactrianMult m (t * s)

-- | Multiplicative proposal with kernel similar to the silhouette of a Bactrian
-- camel. See 'slideBactrian'.
scaleBactrian ::
  SpikeParameter ->
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal Double
scaleBactrian m s = createProposal description (bactrianMultSimple m s) (PDimension 1)
  where
    description = PDescription $ "Scale Bactrian; spike: " ++ show m <> ", sd: " <> show s
