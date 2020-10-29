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
  -- | Forward operator.
  --
  -- For example, for a multiplicative proposal on one variable the forward
  -- operator is @(*)@, so that @x * u = y@.
  (a -> Double -> a) ->
  -- | Inverse operator.
  --
  -- For example, 'recip' for a multiplicative proposal on one variable, since
  -- @y * (recip u) = x * u * (recip u) = x@.
  --
  -- Required for biased proposals.
  Maybe (Double -> Double) ->
  -- | Function to compute the absolute value of the determinant of the Jacobian
  -- matrix. For example, for a multiplicative proposal on one variable, we have
  --
  -- @
  -- detJacobian _ u = Exp $ log $ recip u
  -- @
  --
  -- That is, the determinant of the Jacobian matrix of multiplication is just
  -- the reciprocal value of @u@ (with conversion to log domain).
  --
  -- Required for proposals for which absolute value of the determinant of the
  -- Jacobian differs from 1.0.
  --
  -- Conversion to log domain is necessary, because some determinants of
  -- Jacobians are very small (or large).
  Maybe (a -> Double -> Log Double) ->
  ProposalSimple a
genericContinuous d f mInv mJac x g = do
  u <- genContVar d g
  let r = case mInv of
        Nothing -> 1.0
        Just fInv ->
          let qXY = Exp $ logDensity d u
              qYX = Exp $ logDensity d (fInv u)
           in qYX / qXY
      j = case mJac of
        Nothing -> 1.0
        Just fJac -> fJac x u
  return (x `f` u, r*j)
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
  u <- genDiscreteVar d g
  let r = case mfInv of
        Nothing -> 1.0
        Just fInv ->
          let qXY = Exp $ logProbability d u
              qYX = Exp $ logProbability d (fInv u)
           in qYX / qXY
  return (x `f` u, r)
{-# INLINEABLE genericDiscrete #-}
