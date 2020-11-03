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

    -- * Height labels
    HasHeight (..),
    HeightLabel (..),
    heightL,
    toHeightTree,
    fromHeightTree,
  )
where

import Control.Comonad
import Control.Lens
import Data.Aeson
import Data.Aeson.TH
import Data.Bifunctor
import ELynx.Tree

-- | Should the stem be handled.
--
-- For example, should the stem be considered when calculating the Jacobian
-- during execution of a proposal, or the branch-wise prior?
data HandleStem = WithStem | WithoutStem

-- | Class of types with information about height.
class HasHeight a where
  getHeight :: a -> Length
  setHeight :: Length -> a -> a
  modHeight :: (Length -> Length) -> a -> a

-- | A node label with a height.
--
-- The node height is often used, but height calculation is costly. Direct storage
-- of the node height together with the node label saves time.
newtype HeightLabel a = HeightLabel {fromHeightLabel :: (Length, a)}
  deriving (Show, Eq)

$(deriveJSON defaultOptions ''HeightLabel)

instance HasHeight (HeightLabel a) where
  getHeight = fst . fromHeightLabel
  setHeight x (HeightLabel (_, lb)) = HeightLabel (x, lb)
  modHeight f (HeightLabel (x, lb)) = HeightLabel (f x, lb)

-- | Height of label.
heightL :: HasHeight a => Lens' a Length
heightL = lens getHeight (flip setHeight)

-- | (Re)calculate node heights for a given tree.
toHeightTree :: Tree Length a -> Tree Length (HeightLabel a)
toHeightTree = extend (\t -> HeightLabel (rootHeight t, label t))

-- | Remove information about height from label.
fromHeightTree :: Tree e (HeightLabel a) -> Tree e a
fromHeightTree = second (snd . fromHeightLabel)
