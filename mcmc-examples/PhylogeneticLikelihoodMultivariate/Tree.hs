{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Tree
-- Description :  Markov chain Monte Carlo sampling on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri Jul  3 09:05:09 2020.
--
-- __The import of this module alone should cover most use cases.__
module Tree
  ( T,
    toD,
    fromD,
    getLens,
    nNewick,
    someNewick,
  )
where

import Algebra.Graph.Label
import Algebra.Graph.Labelled.AdjacencyMap
import Data.Aeson
import Data.ByteString.Lazy (ByteString)
import Data.Maybe
import qualified Data.Set as S
import Data.Traversable
import Data.Tree
import Lens.Micro.Platform
import qualified Newick as N

-- | Rooted (directed) tree with branch labels of type @e@ and node labels of
-- type @a@.
type T e a = AdjacencyMap (Distance e) a

instance ToJSON a => ToJSON (Distance a) where
  toJSON = toJSON . fromD
  toEncoding = toEncoding . fromD

instance (FromJSON a, Num a, Ord a) => FromJSON (Distance a) where
  parseJSON d = toD <$> parseJSON d

instance (ToJSON e, ToJSONKey a) => ToJSON (AdjacencyMap e a)

instance (FromJSON e, FromJSONKey a, Ord a) => FromJSON (AdjacencyMap e a)

-- | Extract the number from a distance. Assume that the distance is finite.
fromD :: Distance a -> a
fromD = fromMaybe (error "fromD: Negative distance.") . getFinite . getDistance

-- | Convert a number to a distance. Assume that the given number is positive
-- and finite.
toD :: (Num a, Ord a) => a -> Distance a
toD = distance . fromMaybe (error "toD: Negative number.") . finite

-- | Branch accessor functions (lenses).
getLens :: (Num e, Ord e, Ord a) => a -> a -> Lens' (T e a) e
getLens x y = lens (g x y) (s x y)
  where
    g v w = fromD . edgeLabel v w
    s n m gr e = replaceEdge (toD e) n m gr

label :: Traversable t => t (e, a) -> t (e, (Int, a))
label = snd . mapAccumL (\i (b, x) -> (i + 1, (b, (i, x)))) (1 :: Int)

-- XXX: Removing the origin after insertion is a hack. However, a shared root
-- branch length of zero leads to singular covariance matrix.
fromTree :: (Num e, Ord e, Ord a) => a -> Tree (e, a) -> T e (Int, a)
fromTree o t = removeVertex (0, o) $ (edges . S.toList . go (0, o)) (label t)
  where
    -- go parent tree
    go p (Node (b, x) []) = S.singleton (toD b, p, x)
    go p (Node (b, x) f) = S.insert (toD b, p, x) $ S.unions $ map (go x) f

prune :: (Eq e, Num e) => Tree (e, a) -> Tree (e, a)
prune t@(Node (b, _) [c])
  | b == 0 = prune c
  | otherwise = t
prune t = t

-- | Parse @n@ Newick trees. To ensure unique node labels, 'Int' labels are
-- additionally assigned to all nodes.
nNewick :: Int -> FilePath -> IO [T Double (Int, ByteString)]
nNewick n fn = do
  ts <- N.nNewick n fn
  return $ map (fromTree "origin" . prune) ts

-- | Parse one or more Newick trees until end of file. To ensure unique node
-- labels, 'Int' labels are additionally assigned to all nodes.
someNewick :: FilePath -> IO [T Double (Int, ByteString)]
someNewick fn = do
  ts <- N.someNewick fn
  return $ map (fromTree "origin" . prune) ts
