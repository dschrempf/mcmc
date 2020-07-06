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

import Lens.Micro
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Statistics.Distribution.Normal
import Statistics.Distribution.Uniform

-- The actual proposal with tuning parameter.
slideSimple :: Lens' a Double -> Double -> Double -> Double -> ProposalSimple a
slideSimple l m s t = proposalGenericContinuous l (normalDistr m (s * t)) (+) (-)

-- | Additive proposal with normally distributed density.
slide ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal a
slide n w l m s t = Proposal n w (slideSimple l m s 1.0) tnr
  where
    tnr = if t then Just $ tuner (slideSimple l m s) else Nothing

-- The actual proposal with tuning parameter.
slideSymmetricSimple :: Lens' a Double -> Double -> Double -> ProposalSimple a
slideSymmetricSimple l s t = proposalSymmetricGenericContinuous l (normalDistr 0.0 (s * t)) (+)

-- | Additive proposal with normally distributed density with mean zero. This proposal
-- is very fast, because the Metropolis-Hastings ratio does not include
-- calculation of the forwards and backwards densities.
slideSymmetric ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal a
slideSymmetric n w l s t = Proposal n w (slideSymmetricSimple l s 1.0) tnr
  where
    tnr = if t then Just $ tuner (slideSymmetricSimple l s) else Nothing

-- The actual proposal with tuning parameter.
slideUniformSimple :: Lens' a Double -> Double -> Double -> ProposalSimple a
slideUniformSimple l d t =
  proposalSymmetricGenericContinuous l (uniformDistr (- t * d) (t * d)) (+)

-- | Additive proposal with uniformly distributed density. This proposal is very fast,
-- because the Metropolis-Hastings ratio does not include calculation of the
-- forwards and backwards densities.
slideUniform ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Delta.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal a
slideUniform n w l d t = Proposal n w (slideUniformSimple l d 1.0) tnr
  where
    tnr = if t then Just $ tuner (slideUniformSimple l d) else Nothing
