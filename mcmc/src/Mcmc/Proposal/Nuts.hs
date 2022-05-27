-- |
-- Module      :  Mcmc.Proposal.Nuts
-- Description :  No-U-Turn sampler (NUTS)
-- Copyright   :  (c) 2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri May 27 09:58:23 2022.
--
-- See Matthew D. Hoffman, Andrew Gelman (2014) The No-U-Turn Sampler:
-- Adaptively Setting Path Lengths in Hamiltonian Monte Carlo, Journal of
-- Machine Learning Research.
module Mcmc.Proposal.Nuts
  (
  )
where

import Mcmc.Proposal.Hamiltonian

-- Internal; Slice variable u.
type SliceVariable = Double

-- Internal; Forward is True.
type Direction = Bool

-- Internal; Doubling step number.
type DoublingStep = Int

-- Internal; The number of leapfrog steps within the slice ('n' in Algorithm 3).
type NStepsOK = Int

-- Internal; Stop iteration.
type Stop = Bool

-- Internal; Well, that's fun, isn't it? Have a look at Algorithm 3.
-- cited above.
type BuildTreeReturnType = (Positions, Momenta, Positions, Momenta, Positions, NStepsOK, Stop)

-- NOTE: There are several problems with the NUTS proposal.
--
-- 1. The NUTS proposal is not a real proposal because the new state is to be
-- accepted with probability 1.0. Further, the prior and likelihood are already
-- calculated during the proposal step, and so, should not be recalculated in
-- 'mhgPropose'. A fix requires a (substantial) change of the 'Proposal' data
-- type. All other proposals would also be slower because 'mhgPropose' has to
-- check if the prior and likelihood have been calculated during the proposal
-- step and are already available.
--
-- 2. I do not exactly know how to treat non-unit Jabobians (see
-- 'liftProposalWith'). I think I would need to embed these Jacobians into the
-- posterior calculation of the NUTS proposal. I should also check if this
-- affects the HMC proposal (the prime example is the "Pair" example).

-- buildTreeWith ::
--   Gradient Positions ->
--   Maybe (Validate Positions) ->
--   HMassesInv ->
--   --
--   Positions ->
--   Momenta ->
--   SliceVariable ->
--   Direction ->
--   DoublingStep ->
--   LeapfrogScalingFactor ->
--   BuildTreeReturnType
-- buildTreeWith grad mValidate hMassesInv x p u v j e
--   | j == 0 = case leapfrog grad mValidate hMassesInv 1 e x p of
--       Nothing -> undefined
--       Just _ -> undeefined
--   | otherwise = undefined
