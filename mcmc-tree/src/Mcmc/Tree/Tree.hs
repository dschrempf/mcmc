{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE RankNTypes #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Mcmc.Tree.Tree
-- Description :  Markov chain Monte Carlo sampling on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri Jul  3 09:05:09 2020.
module Mcmc.Tree.Tree
  ( Tree (..),
    fromD,
    toD,
    getLens,
  )
where

import Algebra.Graph.Label
import Algebra.Graph.Labelled.AdjacencyMap
import Data.Aeson
import Data.Maybe
import GHC.Generics
import Lens.Micro

-- | Rooted (directed) tree with node labels of type @a@ and branch labels of
-- type @e@.
newtype Tree e a = Tree {fromTree :: AdjacencyMap (Distance e) a}
  deriving (Generic)

instance ToJSON a => ToJSON (Distance a) where
  toJSON = toJSON . fromD
  toEncoding = toEncoding . fromD

instance (FromJSON a, Num a, Ord a) => FromJSON (Distance a) where
  parseJSON d = toD <$> parseJSON d

instance (ToJSON e, ToJSONKey a) => ToJSON (AdjacencyMap e a)

instance (FromJSON e, FromJSONKey a, Ord a) => FromJSON (AdjacencyMap e a)

instance (ToJSON e, ToJSONKey a) => ToJSON (Tree e a)

instance (FromJSON e, Num e, Ord e, FromJSONKey a, Ord a) => FromJSON (Tree e a)

-- Extract the number from a distance. Assume that the distance is finite.
fromD :: Distance a -> a
fromD = fromMaybe (error "fromD: Negative distance.") . getFinite . getDistance

-- Convert a number to a distance. Assume that the given number is positive
-- and finite.
toD :: (Num a, Ord a) => a -> Distance a
toD = distance . fromMaybe (error "toD: Negative number.") . finite

-- | Branch accessor functions (lenses).
getLens :: (Num e, Ord e, Ord a) => a -> a -> Lens' (Tree e a) e
getLens x y = lens (g x y) (s x y)
  where
    g v w = fromD . edgeLabel v w . fromTree
    s n m gr e = Tree $ replaceEdge (toD e) n m (fromTree gr)
