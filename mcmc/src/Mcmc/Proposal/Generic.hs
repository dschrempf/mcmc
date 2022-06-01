-- |
-- Module      :  Mcmc.Proposal.Generic
-- Description :  Generic interface for creating proposals
-- Copyright   :  (c) Dominik Schrempf 2021
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

-- | Generic function to create proposals for continuous parameters (e.g.,
-- 'Double').
--
-- The procedure is as follows: Let \(\mathbb{X}\) be the state space and \(x\)
-- be the current state.
--
-- 1. Let \(D\) be a continuous probability distribution on \(\mathbb{D}\);
--    sample an auxiliary variable \(epsilon \sim D\).
--
-- 2. Suppose \(\odot : \mathbb{X} \times \mathbb{D} \to \mathbb{X}\). Propose a
--    new state \(x' = x \odot \epsilon\).
--
-- 3. If the proposal is unbiased, the Metropolis-Hastings-Green ratio can
--    directly be calculated using the posterior function.
--
-- 4. However, if the proposal is biased: Suppose \(g : \mathbb{D} \to
--    \mathbb{D}\) inverses the auxiliary variable \(\epsilon\) such that \(x =
--    x' \odot g(\epsilon)\). Calculate the Metropolis-Hastings-Green ratio
--    using the posterior function, \(g\), \(D\), \(\epsilon\), and possibly a
--    Jacobian function.
genericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Probability distribution
  d ->
  -- | Forward operator \(\odot\).
  --
  -- For example, for a multiplicative proposal on one variable the forward
  -- operator is @(*)@, so that @x * u = y@.
  (a -> Double -> a) ->
  -- | Inverse operator \(g\) of the auxiliary variable.
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
  -- Required for proposals for which the absolute value of the determinant of
  -- the Jacobian differs from 1.0.
  --
  -- Conversion to log domain is necessary, because some determinants of
  -- Jacobians are very small (or large).
  Maybe (a -> Double -> Jacobian) ->
  Propose a
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
  pure $ Suggest (x `f` u) r j
{-# INLINEABLE genericContinuous #-}

-- | Generic function to create proposals for discrete parameters (e.g., 'Int').
--
-- See 'genericContinuous'.
genericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Probability distribution.
  d ->
  -- | Forward operator.
  --
  -- For example, (+), so that x + dx = x'.
  (a -> Int -> a) ->
  -- | Inverse operator \(g\) of the auxiliary variable.
  --
  -- For example, 'negate', so that x' + (negate dx) = x.
  --
  -- Only required for biased proposals.
  Maybe (Int -> Int) ->
  Propose a
genericDiscrete d f mfInv x g = do
  u <- genDiscreteVar d g
  let r = case mfInv of
        Nothing -> 1.0
        Just fInv ->
          let qXY = Exp $ logProbability d u
              qYX = Exp $ logProbability d (fInv u)
           in qYX / qXY
  pure $ Suggest (x `f` u) r 1.0
{-# INLINEABLE genericDiscrete #-}
