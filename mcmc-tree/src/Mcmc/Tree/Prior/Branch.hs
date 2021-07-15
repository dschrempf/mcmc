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
  )
where

import Data.List
import ELynx.Tree
import Mcmc.Prior
import Mcmc.Tree.Types

-- | Branch wise prior with given prior function.
branchesWith :: RealFloat a => HandleStem -> PriorFunctionG e a -> PriorFunctionG (Tree e b) a
branchesWith WithStem f (Node br _ ts) = foldl' (*) (f br) $ map (branchesWith WithStem f) ts
branchesWith WithoutStem f (Node _ _ ts) = foldl1' (*) $ map (branchesWith WithStem f) ts
{-# INLINE branchesWith #-}
{-# SPECIALIZE branchesWith :: HandleStem -> PriorFunction e -> PriorFunction (Tree e b) #-}

-- -- NOTE: Somehow I never used the parallelized version, because computing
-- -- priors on branches is just too fast. I leave the implementation for
-- -- further reference.

-- -- | See 'branchesWith'.
-- --
-- -- Evaluate the sub trees up to given layer in parallel. Useful if tree is
-- -- large, or if calculation of the branch prior function is costly.
-- parBranchesWith :: Int -> HandleStem -> PriorFunction e -> PriorFunction (Tree e a)
-- parBranchesWith n WithStem f t = parBranchFoldMap n f (*) t
-- parBranchesWith n WithoutStem f t
--   | n > 1 = foldl1' (*) (map (parBranchFoldMap (n -1) f (*)) ts `using` parList rdeepseq)
--   | otherwise = foldl1' (*) $ map (parBranchFoldMap 0 f (*)) $ forest t
--   where
--     ts = forest t
