{-# LANGUAGE DerivingVia #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Oct 27 19:14:20 2020.
--
-- Type synonyms to improve code readability.

-- |
-- Module      :  Mcmc.Tree.Types
-- Description :  Different tree types
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
module Mcmc.Tree.Types
  ( -- ** Stem
    HandleStem (..),

    -- ** Nodes
    Path,
    HandleNode,
    allNodes,
    withoutRootNode,

    -- ** Ultrametric trees
    toHeightTreeUltrametric,
    fromHeightTree,
  )
where

import ELynx.Tree

-- | Should the stem be handled, when traversing branches of a tree?
data HandleStem = WithStem | WithoutStem

-- | Which nodes should be handled, when traversing a tree?
--
-- Useful when creating proposals on trees.
--
-- For an example, see 'withoutRootNode'.
type HandleNode = Path -> Bool

-- | Handle all nodes.
--
-- In particular:
--
-- - Include the stem, if handling branches.
--
-- - Include the root label, if handling node labels.
allNodes :: HandleNode
allNodes = const True

-- | Exclude the root label.
--
-- @
-- withoutRoot = (>0)
-- @
withoutRootNode :: HandleNode
withoutRootNode = not . null

-- | Calculate node heights for a given tree.
--
-- The __node labels__ and the __stem length__ are __removed__.
--
-- This function is expensive and has not been optimized yet. The run time is
-- @O(n^2)@ where @n@ is the number of inner nodes.
--
-- Return 'Left' if:
--
-- - The tree is not ultrametric. The height of leaves is set to zero. If the
--   tree is not ultrametric, the node heights are not defined and the height
--   tree has to be instantiated manually.
toHeightTreeUltrametric ::
  (HasLength a, Fractional c, Ord c, Show c) =>
  Tree a b ->
  Either String (Tree () c)
-- A leaf.
toHeightTreeUltrametric t
  | ultrametric t = Right $ toHeightTreeUltrametric' t
  | otherwise = Left "toHeightTreeUltrametric: Tree is not ultrametric."

-- Assume the tree is ultrametric.
toHeightTreeUltrametric' :: (HasLength a, Fractional c, Ord c, Show c) => Tree a b -> Tree () c
toHeightTreeUltrametric' t@(Node _ _ ts) =
  Node () (assertNonNegative "toHeightTreeUltrametric'" $ (realToFrac . rootHeight) t) $ map toHeightTreeUltrametric' ts

-- | Remove information about node height from node label.
fromHeightTree :: (Ord a, Num a, Show a) => Tree () a -> Tree a ()
fromHeightTree t = go (label t) t
  where
    go hParent (Node () hNode ts) =
      Node (assertNonNegative "fromHeightTree" (hParent - hNode)) () $ map (go hNode) ts

assertNonNegative :: (Ord a, Num a, Show a) => String -> a -> a
assertNonNegative n val
  | val < 0 = error $ n <> ": Negative value: " <> show val <> "."
  | otherwise = val
