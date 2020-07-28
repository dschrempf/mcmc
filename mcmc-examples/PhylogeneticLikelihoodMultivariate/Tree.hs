{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
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
module Tree
  ( oneTree,
    someTrees,
  )
where

import Codec.Compression.GZip
import Control.Lens
import Data.Aeson
import Data.Bifunctor
import Data.List
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

parseFileWith :: (ShowErrorComponent e) => String -> Parsec e ByteString a -> FilePath -> IO a
parseFileWith s p f = do
  l <- if "gz" `isSuffixOf` f
    then decompress <$> L.readFile f
    else L.readFile f
  return $ case parse p s l of
    Left err -> error $ errorBundlePretty err
    Right val -> val

-- | Parse first Newick tree in file.
oneTree :: FilePath -> IO (Tree Double ByteString)
oneTree p = do
  pt <- parseFileWith "newick" (newick Standard) p
  return $ first fromLength $ either error id $ phyloToLengthTree pt

-- | Parse one or more Newick trees until end of file.
someTrees :: FilePath -> IO [Tree Double ByteString]
someTrees p = do
  pts <- parseFileWith "someNewick" (someNewick Standard) p
  return $ map (first fromLength . either error id . phyloToLengthTree) pts
