{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Tree.Height
-- Description :  Special tree object for time trees
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
  )
where

import Data.Aeson
import Data.Aeson.TH

-- | Class of types with information about height.
class HasHeight a where
  getHeight :: a -> Double
  setHeight :: Double -> a -> a

-- | Change the height.
applyHeight :: HasHeight a => (Double -> Double) -> a -> a
applyHeight f l = setHeight (f $ getHeight l) l

-- | A node label with a height.
--
-- The node height is often used, but height calculation is costly. Direct storage
-- of the node height together with the node label saves time.
newtype HeightLabel a = HeightLabel {fromHeightLabel :: (Double, a)}
  deriving (Show, Eq)

$(deriveJSON defaultOptions ''HeightLabel)

instance HasHeight (HeightLabel a) where
  getHeight = fst . fromHeightLabel
  setHeight x (HeightLabel (_, lb)) = HeightLabel (x, lb)
