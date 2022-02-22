{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal.Slide
-- Description :  Additive proposals
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May  6 10:59:13 2020.
module Mcmc.Proposal.Slide
  ( slide,
    slideSymmetric,
    slideUniformSymmetric,
    slideContrarily,
  )
where

import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Statistics.Types
import Statistics.Distribution.Normal
import Statistics.Distribution.Uniform

-- The actual proposal with tuning parameter.
slideSimple :: Mean Double -> StandardDeviation Double -> TuningParameter -> ProposalSimple Double
slideSimple m s t =
  genericContinuous (normalDistr m (s * t)) (+) (Just negate) Nothing

-- | Additive proposal.
--
-- A normal distribution is used to sample the addend.
slide ::
  Mean Double ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal Double
slide m s = createProposal description (slideSimple m s) (PDimension 1)
  where
    description = PDescription $ "Slide; mean: " ++ show m ++ ", sd: " ++ show s

-- The actual proposal with tuning parameter.
slideSymmetricSimple :: StandardDeviation Double -> TuningParameter -> ProposalSimple Double
slideSymmetricSimple s t =
  genericContinuous (normalDistr 0.0 (s * t)) (+) Nothing Nothing

-- | See 'slide'.
--
-- Use a normal distribution with mean zero. This proposal is fast, because the
-- Metropolis-Hastings-Green ratio does not include calculation of the forwards
-- and backwards kernels.
slideSymmetric ::
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal Double
slideSymmetric s = createProposal description (slideSymmetricSimple s) (PDimension 1)
  where
    description = PDescription $ "Slide symmetric; sd: " ++ show s

-- The actual proposal with tuning parameter.
slideUniformSimple :: Size -> TuningParameter -> ProposalSimple Double
slideUniformSimple d t =
  genericContinuous (uniformDistr (-t * d) (t * d)) (+) Nothing Nothing

-- | See 'slide'.
--
-- Use a uniformly distributed kernel with mean zero. This proposal is fast,
-- because the Metropolis-Hastings-Green ratio does not include calculation of
-- the forwards and backwards kernels.
slideUniformSymmetric ::
  Size ->
  PName ->
  PWeight ->
  Tune ->
  Proposal Double
slideUniformSymmetric d = createProposal description (slideUniformSimple d) (PDimension 1)
  where
    description = PDescription $ "Slide uniform symmetric; delta: " ++ show d

contra :: (Double, Double) -> Double -> (Double, Double)
contra (x, y) u = (x + u, y - u)

slideContrarilySimple ::
  Mean Double ->
  StandardDeviation Double ->
  TuningParameter ->
  ProposalSimple (Double, Double)
slideContrarilySimple m s t =
  genericContinuous (normalDistr m (s * t)) contra (Just negate) Nothing

-- | See 'slide'.
--
-- Use a normally distributed kernel.
--
-- The two values are slid contrarily so that their sum stays constant. Contrary
-- proposals are useful when parameters are confounded.
slideContrarily ::
  Mean Double ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Double, Double)
slideContrarily m s = createProposal description (slideContrarilySimple m s) (PDimension 2)
  where
    description = PDescription $ "Slide contrarily; mean: " ++ show m ++ ", sd: " ++ show s
