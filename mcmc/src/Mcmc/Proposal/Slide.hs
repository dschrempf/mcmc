{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal.Slide
-- Description :  Normally distributed proposal
-- Copyright   :  (c) Dominik Schrempf 2020
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
import Statistics.Distribution.Normal
import Statistics.Distribution.Uniform

-- The actual proposal with tuning parameter.
slideSimple :: Double -> Double -> Double -> ProposalSimple Double
slideSimple m s t =
  genericContinuous (normalDistr m (s * t)) (+) (Just negate) Nothing

-- | Additive proposal with normally distributed kernel.
slide ::
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal Double
slide m s = createProposal description (slideSimple m s)
  where
    description = PDescription $ "Slide; mean: " ++ show m ++ ", sd: " ++ show s

-- The actual proposal with tuning parameter.
slideSymmetricSimple :: Double -> Double -> ProposalSimple Double
slideSymmetricSimple s t =
  genericContinuous (normalDistr 0.0 (s * t)) (+) Nothing Nothing

-- | Additive proposal with normally distributed kernel with mean zero. This
-- proposal is very fast, because the Metropolis-Hastings ratio does not include
-- calculation of the forwards and backwards kernels.
slideSymmetric ::
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal Double
slideSymmetric s = createProposal description (slideSymmetricSimple s)
  where
    description = PDescription $ "Slide symmetric; sd: " ++ show s

-- The actual proposal with tuning parameter.
slideUniformSimple :: Double -> Double -> ProposalSimple Double
slideUniformSimple d t =
  genericContinuous (uniformDistr (- t * d) (t * d)) (+) Nothing Nothing

-- | Additive proposal with uniformly distributed kernel with mean zero. This
-- proposal is very fast, because the Metropolis-Hastings ratio does not include
-- calculation of the forwards and backwards kernels.
slideUniformSymmetric ::
  -- | Delta.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal Double
slideUniformSymmetric d = createProposal description (slideUniformSimple d)
  where
    description = PDescription $ "Slide uniform symmetric; delta: " ++ show d

contra :: (Double, Double) -> Double -> (Double, Double)
contra (x, y) u = (x + u, y - u)

slideContrarilySimple :: Double -> Double -> Double -> ProposalSimple (Double, Double)
slideContrarilySimple m s t =
  genericContinuous (normalDistr m (s * t)) contra (Just negate) Nothing

-- | Additive proposal with normally distributed kernel.
--
-- The two values are slid contrarily so that their sum stays constant. Contrary
-- proposals are useful when parameters are confounded.
slideContrarily ::
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Double, Double)
slideContrarily m s = createProposal description (slideContrarilySimple m s)
  where
    description = PDescription $ "Slide contrarily; mean: " ++ show m ++ ", sd: " ++ show s
