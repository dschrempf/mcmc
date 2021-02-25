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
branchesWith :: HandleStem -> PriorFunction Double -> PriorFunction (Tree Length a)
branchesWith WithStem f = product . map (f . fromLength) . branches
branchesWith WithoutStem f = product . map (f . fromLength) . tail . branches

-- | See 'branchesWith'.
--
-- Evaluate the sub trees up to given layer in parallel. Useful if tree is
-- large, or if the branch prior distribution takes time to evaluate.
parBranchesWith :: Int -> HandleStem -> PriorFunction Double -> PriorFunction (Tree Length a)
parBranchesWith n WithStem f = parBranchFoldMap n (f . fromLength) (*)
parBranchesWith n WithoutStem f = product . map (parBranchFoldMap (n -1) (f . fromLength) (*)) . forest
