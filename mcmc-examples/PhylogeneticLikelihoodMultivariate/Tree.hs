{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
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
  ( oneTree,
    someTrees,
  )
where

import Codec.Compression.GZip
import Control.Lens
-- import Control.Lens.Indexed
-- import Control.Lens.At
import Data.ByteString.Lazy (ByteString)
import qualified Data.ByteString.Lazy as L
import ELynx.Data.Tree
import ELynx.Import.Tree.Newick
import Text.Megaparsec

instance FunctorWithIndex [Int] (Tree e) where
  imap f (Node br lb ts) = Node br (f [] lb) $ imap (\i -> imap (f . (:) i)) ts
  {-# INLINE imap #-}

instance FoldableWithIndex [Int] (Tree e) where
  ifoldMap f (Node _ lb ts) = f [] lb `mappend` ifoldMap (\i -> ifoldMap (f . (:) i)) ts
  {-# INLINE ifoldMap #-}

instance TraversableWithIndex [Int] (Tree e) where
  itraverse f (Node br lb ts) = Node br <$> f [] lb <*> itraverse (\i -> itraverse (f . (:) i)) ts
  {-# INLINE itraverse #-}

type instance IxValue (Tree e a) = (e, a)

type instance Index (Tree e a) = [Int]

instance Ixed (Tree e a) where
  ix xs0 f = go xs0
    where
      go [] (Node br lb ts) = f (br, lb) <&> \(br', lb') -> Node br' lb' ts
      go (i : is) t@(Node br lb ts)
        | i < 0 = pure t
        | otherwise = Node br lb <$> ix i (go is) ts
  {-# INLINE ix #-}

parseFileWith :: (ShowErrorComponent e) => String -> Parsec e ByteString a -> FilePath -> IO a
parseFileWith s p f = do
  l <- decompress <$> L.readFile f
  return $ case parse p s l of
    Left err -> error $ errorBundlePretty err
    Right val -> val

-- | Parse first Newick tree in file.
oneTree :: FilePath -> IO (Tree Length ByteString)
oneTree p = do
  pt <- parseFileWith "newick" (newick Standard) p
  return $ either error id $ phyloToLengthTree pt

-- | Parse one or more Newick trees until end of file.
someTrees :: FilePath -> IO [Tree Length ByteString]
someTrees p = do
  pts <- parseFileWith "someNewick" (someNewick Standard) p
  return $ map (either error id . phyloToLengthTree) pts
