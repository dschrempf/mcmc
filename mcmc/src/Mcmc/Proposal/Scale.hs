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
    scaleContrarily,
  )
where

import Numeric.Log
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Statistics.Distribution.Gamma

-- The actual proposal with tuning parameter. The tuning parameter does not
-- change the mean.
scaleSimple :: Double -> Double -> Double -> ProposalSimple Double
scaleSimple k th t =
  genericContinuous
    (gammaDistr (k / t) (th * t))
    (*)
    (Just recip)
    (Just jac)
  where
    jac _ = Exp . log . recip

-- | Multiplicative proposal with Gamma distributed kernel.
scale ::
  -- | Shape.
  Double ->
  -- | Scale.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Weight ->
  -- | Enable tuning.
  Tune ->
  Proposal Double
scale k th = createProposal description (scaleSimple k th)
  where
    description = "Scale; shape: " ++ show k ++ ", scale: " ++ show th

-- | Multiplicative proposal with Gamma distributed kernel.
--
-- The scale of the Gamma distributions is set to (shape)^{-1}, so that the mean
-- of the Gamma distribution is 1.0.
scaleUnbiased ::
  -- | Shape.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Weight ->
  -- | Enable tuning.
  Tune ->
  Proposal Double
scaleUnbiased k = createProposal description (scaleSimple k (1 / k))
  where
    description = "Scale unbiased; shape: " ++ show k

scaleContrarilySimple :: Double -> Double -> Double -> ProposalSimple (Double, Double)
scaleContrarilySimple k th t =
  genericContinuous
    (gammaDistr (k / t) (th * t))
    contra
    (Just recip)
    (Just jac)
  where contra (x, y) u = (x * u, y / u)
        jac _ u = Exp $ log $ recip $ u*u

-- -- Determinant of Jacobian matrix.
-- contraJac :: (Double, Double) -> Double
-- contraJac (x, y) = x * y


-- | Multiplicative proposal with Gamma distributed kernel.
--
-- The two values are scaled contrarily so that their product stays constant.
-- Contrary proposals are useful when parameters are confounded.
scaleContrarily ::
  -- | Shape.
  Double ->
  -- | Scale.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Weight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Double, Double)
scaleContrarily k th = createProposal description (scaleContrarilySimple k th)
  where
    description = "Scale contrariliy; shape: " ++ show k ++ ", scale: " ++ show th
