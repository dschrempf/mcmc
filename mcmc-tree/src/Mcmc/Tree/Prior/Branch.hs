-- |
-- Module      :  Mcmc.Tree.Prior.Branch
-- Description :  Relative node order constraints and node calibrations
-- Copyright   :  (c) Dominik Schrempf, 2021
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

import Control.Parallel.Strategies
import ELynx.Tree
import Mcmc.Chain.Chain
import Mcmc.Tree.Types

-- | Branch wise prior with given prior function.
branchesWith :: HandleStem -> PriorFunction e -> PriorFunction (Tree e a)
branchesWith = parBranchesWith 0

-- | See 'branchesWith'.
--
-- Evaluate the sub trees up to given layer in parallel. Useful if tree is
-- large, or if calculation of the branch prior function is costly.
parBranchesWith :: Int -> HandleStem -> PriorFunction e -> PriorFunction (Tree e a)
parBranchesWith n WithStem f t = parBranchFoldMap n f (*) t
parBranchesWith n WithoutStem f t
  | n > 1 = product (map (parBranchFoldMap (n-1) f (*)) ts `using` parList rdeepseq)
  | otherwise = product $ map (parBranchFoldMap 0 f (*)) $ forest t
  where ts = forest t
