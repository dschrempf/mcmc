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
    oneTree,
    someTrees,
  )
where

import Algebra.Graph.Label
import AdjacencyIntMap
import Codec.Compression.GZip
import Control.DeepSeq
import Data.Aeson
import Data.Bifunctor
import Data.ByteString.Lazy (ByteString)
import qualified Data.ByteString.Lazy as L
import Data.IntMap (IntMap)
import qualified Data.IntMap as M
import Data.Maybe
import qualified Data.Set as S
import ELynx.Data.Tree
import ELynx.Import.Tree.Newick
import Lens.Micro.Platform
import Text.Megaparsec

-- | Rooted (directed) tree with branch labels of type @e@ and node labels of
-- type @a@.
type T e = AdjacencyIntMap (Distance e)

instance ToJSON e => ToJSON (Distance e) where
  toJSON = toJSON . fromD
  toEncoding = toEncoding . fromD

instance (FromJSON e, Num e, Ord e) => FromJSON (Distance e) where
  parseJSON d = toD <$> parseJSON d

instance (NFData e) => NFData (Distance e) where
  rnf = rnf . getFinite . getDistance

instance ToJSON e => ToJSON (AdjacencyIntMap e)

instance FromJSON e => FromJSON (AdjacencyIntMap e)

-- | Extract the number from a distance. Assume that the distance is finite.
fromD :: Distance e -> e
fromD = fromMaybe (error "fromD: Negative distance.") . getFinite . getDistance

-- | Convert a number to a distance. Assume that the given number is positive
-- and finite.
toD :: (Num e, Ord e) => e -> Distance e
toD = distance . fromMaybe (error "toD: Negative number.") . finite

-- | Branch accessor functions (lenses).
getLens :: (Num e, Ord e) => Int -> Int -> Lens' (T e) e
getLens x y = lens (g x y) (s x y)
  where
    g v w = fromD . edgeLabel v w
    s n m gr e = replaceEdge (toD e) n m gr

-- XXX: Removing the origin after insertion is a hack. However, a shared root
-- branch length of zero leads to singular covariance matrix.
fromTree :: (Num e, Ord e) => Tree e Int -> T e
fromTree t = removeVertex (-1) $ (edges . S.toList . go (-1)) (identify t)
  where
    -- go parent tree
    go p (Node b x []) = S.singleton (toD b, p, x)
    go p (Node b x ts) = S.insert (toD b, p, x) $ S.unions $ map (go x) ts

parseFileWith :: (ShowErrorComponent e) => String -> Parsec e ByteString a -> FilePath -> IO a
parseFileWith s p f = do
  l <- decompress <$> L.readFile f
  return $ case parse p s l of
    Left err -> error $ errorBundlePretty err
    Right val -> val

-- | Parse @n@ Newick trees. Also succeeds when more trees follow.
oneTree :: FilePath -> IO (T Double, IntMap ByteString)
oneTree p = do
  pt <- parseFileWith "oneNewick" (oneNewick Standard) p
  let lt = first fromLength $ either error id $ phyloToLengthTree $ identify pt
      am = fromTree lt
      im = M.fromList $ zip (leaves lt) (leaves pt)
  return (am, im)

-- | Parse one or more Newick trees until end of file.
someTrees :: FilePath -> IO [T Double]
someTrees p = do
  pts <- parseFileWith "someNewick" (someNewick Standard) p
  let lts = map (first fromLength . either error id . phyloToLengthTree . identify) pts
  return $ map fromTree lts
