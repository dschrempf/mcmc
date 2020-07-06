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

import Lens.Micro
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Statistics.Distribution.Gamma

-- The actual proposal with tuning parameter.
scaleSimple :: Lens' a Double -> Double -> Double -> Double -> ProposalSimple a
scaleSimple l k th t = proposalGenericContinuous l (gammaDistr k (t * th)) (*) (/)

-- | Multiplicative proposal with Gamma distributed density.
scale ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Shape.
  Double ->
  -- | Scale.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal a
scale n w l k th t = Proposal n w (scaleSimple l k th 1.0) tnr
  where
    tnr = if t then Just $ tuner $ scaleSimple l k th else Nothing

-- | Multiplicative proposal with Gamma distributed density. The scale of the Gamma
-- distributions is set to (shape)^{-1}, so that the mean of the Gamma
-- distribution is 1.0.
scaleUnbiased ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Shape.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal a
scaleUnbiased n w l k t = Proposal n w (scaleSimple l k (1 / k) 1.0) tnr
  where
    tnr = if t then Just $ tuner $ scaleSimple l k (1 / k) else Nothing
