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
  ( genericContinuous,
    genericDiscrete,
  )
where

import Mcmc.Proposal
import Numeric.Log
import Statistics.Distribution

-- | Generic function to create proposals for continuous parameters ('Double').
genericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Probability distribution
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = x'.
  (a -> Double -> a) ->
  -- | Inverse operator, e.g., 'negate', so that x' + (negate dx) = x. Only
  -- required for biased proposals.
  Maybe (Double -> Double) ->
  ProposalSimple a
genericContinuous d f mfInv x g = do
  dx <- genContVar d g
  let r = case mfInv of
        Nothing -> 1.0
        Just fInv ->
          let qXY = Exp $ logDensity d dx
              qYX = Exp $ logDensity d (fInv dx)
           in qYX / qXY
  return (x `f` dx, r)
{-# INLINEABLE genericContinuous #-}

-- | Generic function to create proposals for discrete parameters ('Int').
genericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Probability distribution.
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = x'.
  (a -> Int -> a) ->
  -- | Inverse operator, e.g., 'negate', so that x' + (negate dx) = x. Only
  -- required for biased proposals.
  Maybe (Int -> Int) ->
  ProposalSimple a
genericDiscrete d f mfInv x g = do
  dx <- genDiscreteVar d g
  let r = case mfInv of
        Nothing -> 1.0
        Just fInv ->
          let qXY = Exp $ logProbability d dx
              qYX = Exp $ logProbability d (fInv dx)
           in qYX / qXY
  return (x `f` dx, r)
{-# INLINEABLE genericDiscrete #-}
