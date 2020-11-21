{-# LANGUAGE TemplateHaskell #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Mcmc.Tree.Types
-- Description :  Different tree types
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Oct 27 19:14:20 2020.
--
-- Type synonyms to improve code readability.
module Mcmc.Tree.Types
  ( -- * Miscellaneous
    HandleStem (..),

    -- * Length trees
    LengthTree,

    -- * Height trees
    HeightLabel (..),
    nodeHeightL,
    nodeNameL,
    HeightTree,
    toHeightTreeUltrametric,
    fromHeightTree,
  )
where

import Control.Lens
import Data.Aeson
import Data.Aeson.TH
import ELynx.Tree

-- | Should the stem be handled?
--
-- For example, should the stem be considered when calculating the Jacobian
-- during execution of a proposal, or the branch-wise prior?
data HandleStem = WithStem | WithoutStem

-- | Length tree.
--
-- A 'Tree' with unconstrained branch lengths and node names.
type LengthTree = Tree Length Name

-- | A node label storing node height and node name.
data HeightLabel = HeightLabel
  { nodeHeight :: Length,
    nodeName :: Name
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''HeightLabel)

-- | Node height.
nodeHeightL :: Lens' HeightLabel Length
nodeHeightL = lens nodeHeight (\x h -> x {nodeHeight = h})

-- | Node name.
nodeNameL :: Lens' HeightLabel Name
nodeNameL = lens nodeName (\x n -> x {nodeName = n})

-- | Height tree.
--
-- A 'Tree' with constrained branch lengths, and node labels storing node height
-- and node name.
--
-- For example, the height tree representation is required when working with
-- ultrametric trees.
type HeightTree = Tree () HeightLabel

-- | Calculate node heights for a given tree.
--
-- __Assumes the tree is ultrametric__ because the height of leaves is set to
-- zero. If the tree is not ultrametric, the node heights cannot be obtained and
-- the height tree has to be instantiated manually.
--
-- The length of the stem is lost.
--
-- This operation is slow, O(n^2), where n is the number of inner nodes.
toHeightTreeUltrametric :: LengthTree -> HeightTree
toHeightTreeUltrametric t@(Node _ lb ts) =
  Node
    ()
    (HeightLabel (rootHeight t) lb)
    (map toHeightTreeUltrametric ts)

-- | Remove information about node height from node label.
fromHeightTree :: HeightTree -> LengthTree
fromHeightTree t = go (nodeHeight $ label t) t
  where
    go hParent (Node () lb ts) =
      let hNode = nodeHeight lb
          nNode = nodeName lb
       in Node (hParent - hNode) nNode $ map (go hNode) ts

-- -- Use the height of the parent node to set the stem length of the tree.
-- setStemLength :: Length -> HeightTree -> HeightTree
-- setStemLength hParent t = t & branchL .~ dh
--   where
--     hChild = t ^. labelL . measurableL
--     dh = hParent - hChild

-- -- | Recalculate the branch lengths.
-- --
-- -- The stem length is unchanged.
-- --
-- -- This operation is slow.
-- recalculateBranchLengths :: HeightTree -> HeightTree
-- recalculateBranchLengths t =
--   over (forestL . mapped) (setStemLength l . recalculateBranchLengths) t
--   where
--     l = t ^. labelL . measurableL

-- -- | Recalculate the branch lengths up to a given depth.
-- --
-- -- We have:
-- -- @
-- -- recalculateBranchLengthsNLevels 0 t = t
-- -- @
-- --
-- -- The stem length is unchanged.
-- --
-- -- This operation is slow.
-- recalculateBranchLengthsNLevels :: Int -> HeightTree -> HeightTree
-- recalculateBranchLengthsNLevels 0 t = t
-- recalculateBranchLengthsNLevels n t =
--   over (forestL . mapped) (setStemLength l . recalculateBranchLengthsNLevels (n -1)) t
--   where
--     l = t ^. labelL . measurableL
