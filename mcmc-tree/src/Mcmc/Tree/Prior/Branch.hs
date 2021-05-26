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
--
-- Branch wise priors.
module Mcmc.Tree.Prior.Branch
  ( branchesWith,
    parBranchesWith,
  )
where

import ELynx.Tree
import Mcmc.Chain.Chain
import Mcmc.Tree.Types

-- | Branch wise prior with given prior function.
branchesWith :: HandleDepth -> PriorFunction e -> PriorFunction (Tree e a)
branchesWith = parBranchesWith 0

-- | See 'branchesWith'.
--
-- Evaluate the sub trees up to given layer in parallel. Useful if tree is
-- large, or if the branch prior distribution takes time to evaluate.
parBranchesWith :: Int -> HandleDepth -> PriorFunction e -> PriorFunction (Tree e a)
parBranchesWith n hd f = parBranchFoldMapWithDepth n f' (*)
  where f' d br = if hd d then f br else 1.0
