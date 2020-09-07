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

-- | Branch length prior with given distribution.
--
-- Parallel version; evaluate the sub trees of the root node in parallel. Only
-- suitable for balanced trees with thousands of leaves.
parBranchesWith ::
  -- | Branch prior distribution.
  (Double -> Log Double) ->
  Tree Double a ->
  Log Double
parBranchesWith f = parBranchFoldMap 1 f (*)
