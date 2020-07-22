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
    proposalSymmetricGenericContinuous,
    proposalGenericDiscrete,
    proposalSymmetricGenericDiscrete,
  )
where

import Mcmc.Proposal
import Numeric.Log
import Statistics.Distribution
import System.Random.MWC

jumpCont ::
  (ContDistr d, ContGen d) =>
  d ->
  (Double -> Double -> Double) ->
  Double ->
  GenIO ->
  IO Double
jumpCont d f x g = do
  dx <- genContVar d g
  return $ x `f` dx
{-# INLINEABLE jumpCont #-}

kernelCont ::
  (ContDistr d, ContGen d) =>
  d ->
  (Double -> Double -> Double) ->
  Double ->
  Double ->
  Log Double
kernelCont d fInv x y = Exp $ logDensity d (y `fInv` x)
{-# INLINEABLE kernelCont #-}

-- | Generic function to create proposals for continuous parameters ('Double').
proposalGenericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Probability distribution
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Double -> Double -> Double) ->
  -- | Inverse operator, e.g.,(-), so that y - dx = x.
  (Double -> Double -> Double) ->
  ProposalSimple Double
proposalGenericContinuous d f fInv =
  ProposalSimple (jumpCont d f) (Just $ kernelCont d fInv)

-- | Generic function to create symmetric proposals for continuous parameters ('Double').
proposalSymmetricGenericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Probability distribution
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Double -> Double -> Double) ->
  ProposalSimple Double
proposalSymmetricGenericContinuous d f =
  ProposalSimple (jumpCont d f) Nothing

jumpDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  d ->
  (Int -> Int -> Int) ->
  Int ->
  GenIO ->
  IO Int
jumpDiscrete d f x g = do
  dx <- genDiscreteVar d g
  return $ x `f` dx
{-# INLINEABLE jumpDiscrete #-}

kernelDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  d ->
  (Int -> Int -> Int) ->
  Int ->
  Int ->
  Log Double
kernelDiscrete d fInv x y =
  Exp $ logProbability d (y `fInv` x)
{-# INLINEABLE kernelDiscrete #-}

-- | Generic function to create proposals for discrete parameters ('Int').
proposalGenericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Probability distribution.
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Int -> Int -> Int) ->
  -- | Inverse operator, e.g.,(-), so that y - dx = x.
  (Int -> Int -> Int) ->
  ProposalSimple Int
proposalGenericDiscrete fd f fInv =
  ProposalSimple (jumpDiscrete fd f) (Just $ kernelDiscrete fd fInv)

-- | Generic function to create symmetric proposals for discrete parameters ('Int').
proposalSymmetricGenericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Probability distribution.
  d ->
  -- | Forward operator, e.g. (+), so that x + dx = y.
  (Int -> Int -> Int) ->
  ProposalSimple Int
proposalSymmetricGenericDiscrete fd f =
  ProposalSimple (jumpDiscrete fd f) Nothing
