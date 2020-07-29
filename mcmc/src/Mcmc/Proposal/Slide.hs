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
    slideUniform,
  )
where

import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Statistics.Distribution.Normal
import Statistics.Distribution.Uniform

-- The actual proposal with tuning parameter.
slideSimple :: Double -> Double -> Double -> ProposalSimple Double
slideSimple m s t = proposalGenericContinuous (normalDistr m (s * t)) (+) (Just negate)

-- | Additive proposal with normally distributed kernel.
slide ::
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
  Proposal Double
slide n w m s = createProposal n w (slideSimple m s)

-- The actual proposal with tuning parameter.
slideSymmetricSimple :: Double -> Double -> ProposalSimple Double
slideSymmetricSimple s t = proposalGenericContinuous (normalDistr 0.0 (s * t)) (+) Nothing

-- | Additive proposal with normally distributed kernel with mean zero. This
-- proposal is very fast, because the Metropolis-Hastings ratio does not include
-- calculation of the forwards and backwards kernels.
slideSymmetric ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal Double
slideSymmetric n w s = createProposal n w (slideSymmetricSimple s)

-- The actual proposal with tuning parameter.
slideUniformSimple :: Double -> Double -> ProposalSimple Double
slideUniformSimple d t =
  proposalGenericContinuous (uniformDistr (- t * d) (t * d)) (+) Nothing

-- | Additive proposal with uniformly distributed kernel. This proposal is very fast,
-- because the Metropolis-Hastings ratio does not include calculation of the
-- forwards and backwards kernels.
slideUniform ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Delta.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal Double
slideUniform n w d = createProposal n w (slideUniformSimple d)
