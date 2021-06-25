{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingVia #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TypeFamilies #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Mcmc.Tree.Types
-- Description :  Different tree types
-- Copyright   :  (c) Dominik Schrempf, 2021
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
  ( -- ** Stem
    HandleStem (..),

    -- ** Nodes
    Path,
    HandleNode,
    allNodes,
    withoutRootNode,

    -- ** Heights
    Height (fromHeight),
    HasHeight (..),
    toHeight,
    toHeightUnsafe,
    HeightLabel (..),
    nodeHeightL,
    nodeNameL,

    -- ** Height trees
    HeightTree,
    toHeightTreeUltrametric,
    fromHeightTree,
  )
where

import Control.DeepSeq
import Control.Lens
import Data.Aeson
import Data.Aeson.TH
import Data.Monoid
import Data.Vector.Unboxed.Deriving
import ELynx.Tree
import GHC.Generics

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

-- | Non-negative height.
--
-- However, non-negativity is only checked with 'toHeight', and negative values
-- can be obtained using the 'Num' and related instances.
--
-- Safe conversion is roughly 50 percent slower.
newtype Height = Height {fromHeight :: Double}
  deriving (Read, Show, Generic, NFData)
  deriving (Enum, Eq, Floating, Fractional, Num, Ord, Real, RealFloat, RealFrac) via Double
  deriving (Semigroup, Monoid) via Sum Double

derivingUnbox
  "Height"
  [t|Height -> Double|]
  [|fromHeight|]
  [|Height|]

$(deriveJSON defaultOptions ''Height)

-- | If negative, call 'error' with given calling function name.
toHeight :: Double -> Either String Height
toHeight x
  | x < 0 = Left $ "toHeight: Height is negative: " ++ show x ++ "."
  | otherwise = Right $ Height x

-- | Do not check if value is negative.
toHeightUnsafe :: Double -> Height
toHeightUnsafe = Height

-- | A data type with measurable and modifiable values.
class HasHeight a where
  getHeight :: a -> Height
  setHeight :: Height -> a -> a
  modHeight :: (Height -> Height) -> a -> a

-- | A node label storing node height and node name.
data HeightLabel a = HeightLabel
  { nodeHeight :: Height,
    nodeName :: a
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''HeightLabel)

instance HasHeight (HeightLabel a) where
  getHeight = nodeHeight
  setHeight h (HeightLabel _ lb) = HeightLabel h lb
  modHeight f (HeightLabel h lb) = HeightLabel (f h) lb

-- | Node height.
nodeHeightL :: Lens' (HeightLabel a) Height
nodeHeightL = lens nodeHeight (\x h -> x {nodeHeight = h})

-- | Node name.
nodeNameL :: Lens' (HeightLabel a) a
nodeNameL = lens nodeName (\x n -> x {nodeName = n})

-- | Height tree.
--
-- A 'Tree' with constrained branch lengths, and node labels storing node height
-- and node name.
--
-- For example, the height tree representation is required when working with
-- ultrametric trees.
type HeightTree a = Tree () (HeightLabel a)

-- | Calculate node heights for a given tree.
--
-- The __length of the stem is lost__.
--
-- This function is expensive and has not been optimized yet. The run time is
-- @O(n^2)@ where @n@ is the number of inner nodes.
--
-- Return 'Left' if:
--
-- - The tree is not ultrametric. The height of leaves is set to zero. If the
--   tree is not ultrametric, the node heights are not defined and the height
--   tree has to be instantiated manually.
toHeightTreeUltrametric :: Tree Length a -> Either String (HeightTree a)
-- A leaf.
toHeightTreeUltrametric t
  | ultrametric t = Right $ toHeightTreeUltrametric' t
  | otherwise = Left "toHeightTreeUltrametric: Tree is not ultrametric."

-- Assume the tree is ultrametric.
toHeightTreeUltrametric' :: Tree Length a -> HeightTree a
toHeightTreeUltrametric' t@(Node _ lb ts) =
  Node
    ()
    ( HeightLabel
        ( either (error . (<>) "toHeightTreeUltrametric': ") id $
            toHeight $ fromLength $ rootHeight t
        )
        lb
    )
    (map toHeightTreeUltrametric' ts)

-- | Remove information about node height from node label.
fromHeightTree :: HeightTree a -> Tree Length a
fromHeightTree t = go (nodeHeight $ label t) t
  where
    go hParent (Node () lb ts) =
      let hNode = nodeHeight lb
          nNode = nodeName lb
       in Node
            ( either (error . (<>) "fromHeightTree: ") id $
                toLength $ fromHeight $ hParent - hNode
            )
            nNode
            $ map (go hNode) ts
