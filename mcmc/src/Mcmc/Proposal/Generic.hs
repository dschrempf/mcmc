{-# LANGUAGE RankNTypes #-}

-- Technically, only a Getter is needed when calculating the density of the proposal
-- ('densCont', and similar functions). I tried splitting the lens into a getter
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
    proposalSymmetricGenericContinuous,
    proposalGenericDiscrete,
    proposalSymmetricGenericDiscrete,
  )
where

import Lens.Micro
import Mcmc.Proposal
import Numeric.Log
import Statistics.Distribution
import System.Random.MWC

jumpCont ::
  (ContDistr d, ContGen d) =>
  Lens' a Double ->
  d ->
  (Double -> Double -> Double) ->
  a ->
  GenIO ->
  IO a
jumpCont l d f x g = do
  dx <- genContVar d g
  return $ set l ((x ^. l) `f` dx) x
{-# INLINEABLE jumpCont #-}

densCont ::
  (ContDistr d, ContGen d) =>
  Lens' a Double ->
  d ->
  (Double -> Double -> Double) ->
  a ->
  a ->
  Log Double
densCont l d fInv x y = Exp $ logDensity d ((y ^. l) `fInv` (x ^. l))
{-# INLINEABLE densCont #-}

-- | Generic function to create proposals for continuous parameters ('Double').
proposalGenericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Probability distribution
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Double -> Double -> Double) ->
  -- | Inverse operator, e.g.,(-), so that y - dx = x.
  (Double -> Double -> Double) ->
  ProposalSimple a
proposalGenericContinuous l d f fInv =
  ProposalSimple (jumpCont l d f) (Just $ densCont l d fInv)

-- | Generic function to create symmetric proposals for continuous parameters ('Double').
proposalSymmetricGenericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Instruction about which parameter to change.
  Lens' a Double ->
  -- | Probability distribution
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Double -> Double -> Double) ->
  ProposalSimple a
proposalSymmetricGenericContinuous l d f =
  ProposalSimple (jumpCont l d f) Nothing

jumpDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  Lens' a Int ->
  d ->
  (Int -> Int -> Int) ->
  a ->
  GenIO ->
  IO a
jumpDiscrete l d f x g = do
  dx <- genDiscreteVar d g
  return $ set l ((x ^. l) `f` dx) x
{-# INLINEABLE jumpDiscrete #-}

densDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  Lens' a Int ->
  d ->
  (Int -> Int -> Int) ->
  a ->
  a ->
  Log Double
densDiscrete l d fInv x y =
  Exp $ logProbability d ((y ^. l) `fInv` (x ^. l))
{-# INLINEABLE densDiscrete #-}

-- | Generic function to create proposals for discrete parameters ('Int').
proposalGenericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Instruction about which parameter to change.
  Lens' a Int ->
  -- | Probability distribution.
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Int -> Int -> Int) ->
  -- | Inverse operator, e.g.,(-), so that y - dx = x.
  (Int -> Int -> Int) ->
  ProposalSimple a
proposalGenericDiscrete l fd f fInv =
  ProposalSimple (jumpDiscrete l fd f) (Just $ densDiscrete l fd fInv)

-- | Generic function to create symmetric proposals for discrete parameters ('Int').
proposalSymmetricGenericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Instruction about which parameter to change.
  Lens' a Int ->
  -- | Probability distribution.
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Int -> Int -> Int) ->
  ProposalSimple a
proposalSymmetricGenericDiscrete l fd f =
  ProposalSimple (jumpDiscrete l fd f) Nothing
