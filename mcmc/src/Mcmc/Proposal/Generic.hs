{-# LANGUAGE RankNTypes #-}

-- Technically, only a Getter is needed when calculating the kernel of the proposal
-- ('kernelCont', and similar functions). I tried splitting the lens into a getter
-- and a setter. However, speed improvements were marginal, and some times not
-- even measurable. Using a 'Lens'' is just easier, and has no real drawbacks.

-- |
-- Module      :  Mcmc.Proposal.Generic
-- Description :  Generic interface to create proposals
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu May 14 20:26:27 2020.
module Mcmc.Proposal.Generic
  ( proposalGenericContinuous,
    proposalGenericDiscrete,
  )
where

import Mcmc.Proposal
import Numeric.Log
import Statistics.Distribution
import System.Random.MWC

sampleCont ::
  (ContDistr d, ContGen d) =>
  d ->
  (a -> Double -> a) ->
  Maybe (Double -> Double) ->
  a ->
  GenIO ->
  IO (a, Log Double)
sampleCont d f mfInv x g = do
  dx <- genContVar d g
  case mfInv of
    Nothing -> return (x `f` dx, 1.0)
    Just fInv -> do
      let dxInv = fInv dx
          qXY = Exp $ logDensity d dx
          qYX = Exp $ logDensity d dxInv
      return (x `f` dx, qYX / qXY)
{-# INLINEABLE sampleCont #-}

-- | Generic function to create proposals for continuous parameters ('Double').
proposalGenericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Probability distribution
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (a -> Double -> a) ->
  -- | Inverse operator, e.g., 'negate' or 'recip'; only required for biased
  -- proposals.
  Maybe (Double -> Double) ->
  ProposalSimple a
proposalGenericContinuous d f fInv = ProposalSimple (sampleCont d f fInv)

sampleDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  d ->
  (a -> Int -> a) ->
  Maybe (Int -> Int) ->
  a ->
  GenIO ->
  IO (a, Log Double)
sampleDiscrete d f mfInv x g = do
  dx <- genDiscreteVar d g
  case mfInv of
    Nothing -> return (x `f` dx, 1.0)
    Just fInv -> do
      let dxInv = fInv dx
          qXY = Exp $ logProbability d dx
          qYX = Exp $ logProbability d dxInv
      return (x `f` dx, qYX / qXY)
{-# INLINEABLE sampleDiscrete #-}

-- | Generic function to create proposals for discrete parameters ('Int').
proposalGenericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Probability distribution.
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (a -> Int -> a) ->
  -- | Inverse operator, e.g., 'negate'; only required for biased proposals.
  Maybe (Int -> Int) ->
  ProposalSimple a
proposalGenericDiscrete fd f fInv = ProposalSimple (sampleDiscrete fd f fInv)
