{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal.Scale
-- Description :  Multiplicative proposals
-- Copyright   :  (c) Dominik Schrempf 2021
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

import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Statistics.Types
import Numeric.Log
import Statistics.Distribution.Gamma

-- The actual proposal with tuning parameter. The tuning parameter does not
-- change the mean.
scaleSimple :: Shape Double -> Scale Double -> TuningParameter -> ProposalSimple Double
scaleSimple k th t =
  genericContinuous
    (gammaDistr (k / t) (th * t))
    (*)
    (Just recip)
    (Just jac)
  where
    jac _ = Exp . log . recip

-- | Multiplicative proposal.
--
-- The gamma distribution is used to sample the multiplier. Therefore, this and
-- all derived proposals are log-additive in that they do not change the sign of
-- the state. Further, the value zero is never proposed when having a strictly
-- positive value.
--
-- Consider using 'Mcmc.Proposal.Slide.slide' to allow proposition of values
-- having opposite sign.
scale ::
  Shape Double ->
  Scale Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal Double
scale k th = createProposal description (scaleSimple k th) (PDimension 1)
  where
    description = PDescription $ "Scale; shape: " ++ show k ++ ", scale: " ++ show th

-- | See 'scale'.
--
-- The scale of the gamma distribution is set to (shape)^{-1}, so that the mean
-- of the gamma distribution is 1.0.
scaleUnbiased ::
  Shape Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal Double
scaleUnbiased k = createProposal description (scaleSimple k (1 / k)) (PDimension 1)
  where
    description = PDescription $ "Scale unbiased; shape: " ++ show k

scaleContrarilySimple ::
  Shape Double ->
  Scale Double ->
  TuningParameter ->
  ProposalSimple (Double, Double)
scaleContrarilySimple k th t =
  genericContinuous
    (gammaDistr (k / t) (th * t))
    contra
    (Just recip)
    (Just jac)
  where
    contra (x, y) u = (x * u, y / u)
    jac _ u = Exp $ log $ recip $ u * u

-- | See 'scale'.
--
-- The two values are scaled contrarily so that their product stays constant.
-- Contrary proposals are useful when parameters are confounded.
scaleContrarily ::
  Shape Double ->
  Scale Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Double, Double)
scaleContrarily k th = createProposal description (scaleContrarilySimple k th) (PDimension 2)
  where
    description = PDescription $ "Scale contrariliy; shape: " ++ show k ++ ", scale: " ++ show th
