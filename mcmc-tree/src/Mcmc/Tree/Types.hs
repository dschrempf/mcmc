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
  ( SubstitutionTree,
    TimeTree,
    toTimeTree,
    fromTimeTree,
    RateTree,
  )
where

import Control.Comonad
import Data.Aeson
import Data.Bifunctor
import qualified Data.ByteString.Char8 as BS
import ELynx.Tree
import Mcmc.Tree.Height

-- | Substitution tree.
--
-- The branches are measured in number of substitutions.
--
-- The node labels store the node names.
type SubstitutionTree = Tree Double BS.ByteString

-- | Time tree.
--
-- The branches are measured in time and denote durations.
--
-- The node labels store the node ages together with the names.
type TimeTree = Tree Double (HeightLabel BS.ByteString)

-- | Calculate node ages.
toTimeTree :: Tree Double BS.ByteString -> TimeTree
toTimeTree = extend (\t -> HeightLabel (rootHeight t, label t))

-- | Forget node ages.
fromTimeTree :: TimeTree -> Tree Double BS.ByteString
fromTimeTree = second (snd . fromHeightLabel)

-- | Rate tree.
--
-- The branches are measured in relative or absolute rates.
--
-- The node labels store names.
type RateTree = Tree Double BS.ByteString

-- This is pretty lame, but I need JSON instances.

instance ToJSON BS.ByteString where
  toJSON = toJSON . BS.unpack
  toEncoding = toEncoding . BS.unpack

instance FromJSON BS.ByteString where
  parseJSON = fmap BS.pack . parseJSON
