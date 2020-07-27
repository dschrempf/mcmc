{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal.Scale
-- Description :  Scaling proposal with Gamma distribution
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu May 14 21:49:23 2020.
module Mcmc.Proposal.Scale
  ( scale,
    scaleUnbiased,
  )
where

import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Statistics.Distribution.Gamma

-- The actual proposal with tuning parameter. The tuning parameter does not
-- change the mean.
scaleSimple :: Double -> Double -> Double -> ProposalSimple Double
scaleSimple k th t = proposalGenericContinuous (gammaDistr (k / t) (t * th)) (*) (/)

-- | Multiplicative proposal with Gamma distributed kernel.
scale ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Scale.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal Double
scale n w k th t = Proposal n w (scaleSimple k th 1.0) tnr
  where
    tnr = if t then Just $ tuner $ scaleSimple k th else Nothing

-- | Multiplicative proposal with Gamma distributed kernel. The scale of the Gamma
-- distributions is set to (shape)^{-1}, so that the mean of the Gamma
-- distribution is 1.0.
scaleUnbiased ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal Double
scaleUnbiased n w k t = Proposal n w (scaleSimple k (1 / k) 1.0) tnr
  where
    tnr = if t then Just $ tuner $ scaleSimple k (1 / k) else Nothing
