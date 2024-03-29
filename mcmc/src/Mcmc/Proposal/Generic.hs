-- |
-- Module      :  Mcmc.Proposal.Generic
-- Description :  Generic interface for creating proposals
-- Copyright   :  2021 Dominik Schrempf
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

-- | Generic function to create proposals using a continuous auxiliary variable
-- of type 'Double'.
--
-- The procedure is as follows: Let \(\mathbb{X}\) be the state space and \(x\)
-- be the current state.
--
-- 1. Let \(D\) be a continuous probability distribution on \(\mathbb{D}\);
--    sample an auxiliary variable \(u \sim D\).
--
-- 2. Let \(\odot : \mathbb{X} \times \mathbb{D} \to \mathbb{X}\). Propose a
--    new state \(x' = x \odot u\).
--
-- If the proposal is unbiased, the Metropolis-Hastings-Green ratio can directly
-- be calculated using the posterior function.
--
-- However, if the proposal is biased: Let \(g : \mathbb{D} \to \mathbb{D}\);
-- \(g\) inverses the auxiliary variable \(u\) such that \(x = x' \odot g(u)\).
-- Calculate the Metropolis-Hastings-Green ratio using the posterior function,
-- \(g\), \(D\), \(u\), and possibly a Jacobian function.
genericContinuous ::
  (ContDistr d, ContGen d) =>
  -- | Probability distribution
  d ->
  -- | Forward operator \(\odot\).
  --
  -- For example, for a multiplicative proposal on one variable the forward
  -- operator is @(*)@, so that \(x' = x * u\).
  (a -> Double -> a) ->
  -- | Inverse operator \(g\) of the auxiliary variable.
  --
  -- For example, 'recip' for a multiplicative proposal on one variable, since
  -- \(x' * u^{-1} = x * u * u^{-1} = x\).
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
  PFunction a
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
  pure (Propose (x `f` u) r j, Nothing)
{-# INLINEABLE genericContinuous #-}

-- | Generic function to create proposals using a discrete auxiliary variable of
-- type 'Int'.
--
-- See 'genericContinuous'.
genericDiscrete ::
  (DiscreteDistr d, DiscreteGen d) =>
  -- | Probability distribution.
  d ->
  -- | Forward operator.
  --
  -- For example, (+), so that \(x + dx = x'\).
  (a -> Int -> a) ->
  -- | Inverse operator \(g\) of the auxiliary variable.
  --
  -- For example, 'negate', so that \(x' - dx = x + dx - dx = x\).
  --
  -- Only required for biased proposals.
  Maybe (Int -> Int) ->
  PFunction a
genericDiscrete d f mfInv x g = do
  u <- genDiscreteVar d g
  let r = case mfInv of
        Nothing -> 1.0
        Just fInv ->
          let qXY = Exp $ logProbability d u
              qYX = Exp $ logProbability d (fInv u)
           in qYX / qXY
  pure (Propose (x `f` u) r 1.0, Nothing)
{-# INLINEABLE genericDiscrete #-}
