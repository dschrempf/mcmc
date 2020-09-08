-- |
-- Module      :  Mcmc.Tree.Prior.Branch
-- Description :  Relative node order constraints and node calibrations
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jul 27 10:49:11 2020.
module Mcmc.Tree.Prior.Branch
  ( branchesWith,
    parBranchesWith,
  )
where

import ELynx.Tree
import Numeric.Log

-- | Branch length prior with given distribution.
branchesWith ::
  -- | Branch prior distribution.
  (Double -> Log Double) ->
  Tree Double a ->
  Log Double
branchesWith f = product . map f . branches

-- | See 'branchesWith'.
--
-- Evaluate the sub trees up to given layer in Useful if tree is large, or if
-- the branch prior distribution takes time to evaluate.
parBranchesWith ::
  Int ->
  (Double -> Log Double) ->
  Tree Double a ->
  Log Double
parBranchesWith n f = parBranchFoldMap n f (*)
