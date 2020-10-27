-- |
-- Module      :  Mcmc.Tree.Height
-- Description :  Special tree object storing height
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Oct 27 17:30:47 2020.
module Mcmc.Tree.Height
  ( HasHeight (..),
    applyHeight,
    HeightLabel (..),
    HeightTree,
  )
where

import ELynx.Tree

-- | A class for node labels that have an associated height.
class HasHeight a where
  getHeight :: a -> Double
  setHeight :: Double -> a -> a

-- | Change the height.
applyHeight :: HasHeight a => (Double -> Double) -> a -> a
applyHeight f l = setHeight (f $ getHeight l) l

-- | Directly store the node height together with the node label.
--
-- The node height is often used and queried repeatedly, but height calculation
-- is costly. Here, the node height is stored together with the node label.
newtype HeightLabel a = HeightLabel {getLabel :: (Double, a)}
  deriving (Show, Eq)

instance HasHeight (HeightLabel a) where
  getHeight = fst . getLabel
  setHeight x (HeightLabel (_, lb)) = HeightLabel (x, lb)

-- | A tree with stored node height.
type HeightTree a = Tree Double (HeightLabel a)
