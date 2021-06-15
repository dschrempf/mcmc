-- |
-- Module      :  Mcmc.Tree.Mrca
-- Description :  Work with paths to nodes on trees
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Oct 27 12:55:07 2020.
module Mcmc.Tree.Mrca
  ( -- ** Most recent command ancestors
    mrca,
  )
where

import Data.Either
import qualified Data.Set as S
import ELynx.Tree

-- Test if the root node of the given tree is an ancestor of the given leaves.
isAncestor :: Ord a => S.Set a -> Tree e a -> Bool
--                      True if an x of xs is not in the collection of leaves.
--                False      if an x of xs is not in the collection of leaves. -> OK.
isAncestor xs t = not $ any (`S.notMember` lvs) xs
  where
    lvs = S.fromList $ leaves t

-- Test if the root node of the given tree is the MRCA of the given leaves.
isMrca :: Ord a => S.Set a -> Tree e a -> Bool
--                                    True if any daughter forest is an ancestor.
--                               False     if any daughter forest is an ancestor. -> OK.
isMrca xs t = isAncestor xs t && not (any (isAncestor xs) (forest t))

-- | Get the path to the MRCA of the given list of leaves on the given tree.
--
-- Return 'Left' if:
--
-- - The leaves of the tree contain duplicates.
--
-- - The list of leaves contains duplicates.
--
-- - The MRCA cannot be found.
mrca :: (Ord a, Show a) => [a] -> Tree e a -> Either String Path
mrca xs tr
  | duplicateLeaves tr = Left "mrca: Tree contains duplicate leaves."
  | S.size ss /= length xs = Left $ "mrca: List of leaves contains duplicates: " <> show xs <> "."
  | otherwise = tail <$> go 0 tr
  where
    ss = S.fromList xs
    go i t
      | isMrca ss t = Right [i]
      --                                    One path will be (Right p).
      | isAncestor ss t = Right $ i : head (rights [go j t' | (j, t') <- zip [0 ..] (forest t)])
      | otherwise = Left $ "Could not get MRCA for: " <> show xs <> "."
