-- |
-- Module      :  Mcmc.Tree.Mrca
-- Description :  Work with paths to nodes on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Oct 27 12:55:07 2020.
module Mcmc.Tree.Mrca
  ( mrca,
    mrcaUnsafe,
  )
where

import Data.Either
import qualified Data.Set as S
import ELynx.Tree

-- Test if the root node of the given tree is an ancestor of the given leaves.
isAncestor :: Ord a => [a] -> Tree e a -> Bool
--                      True if an x of xs is not in the collection of leaves.
--                False      if an x of xs is not in the collection of leaves. -> OK.
isAncestor xs t = not $ any (`S.notMember` lvs) xs
  where
    lvs = S.fromList $ leaves t

-- Test if the root node of the given tree is the MRCA of the given leaves.
isMrca :: Ord a => [a] -> Tree e a -> Bool
--                                    True if any daughter forest is an ancestor.
--                               False     if any daughter forest is an ancestor. -> OK.
isMrca xs t = isAncestor xs t && not (any (isAncestor xs) (forest t))

-- | Get the path to the MRCA of the given list of node labels on the given
-- tree.
--
-- Assume that the leaves are unique.
mrca :: (Ord a, Show a) => [a] -> Tree e a -> Either String Path
mrca xs = fmap tail . go 0
  where
    go i t
      | isMrca xs t = Right [i]
      --                                    One path will be (Right p).
      | isAncestor xs t = Right $ i : head (rights [go j t' | (j, t') <- zip [0 ..] (forest t)])
      | otherwise = Left $ "Could not get MRCA for: " <> show xs

-- | See 'mrca'.
--
-- Call 'error' if the MRCA is not found on the tree.
mrcaUnsafe :: (Ord a, Show a) => [a] -> Tree e a -> Path
mrcaUnsafe xs = either error id . mrca xs
