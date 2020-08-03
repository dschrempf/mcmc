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

import Data.Attoparsec.Lazy
import Codec.Compression.GZip
import Control.Lens
import Data.Bifunctor
import Data.List
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as BL
import ELynx.Data.Tree
import ELynx.Import.Tree.Newick

instance FunctorWithIndex [Int] (Tree e) where
  imap f (Node br lb ts) = Node br (f [] lb) $ imap (\i -> imap (f . (:) i)) ts
  {-# INLINE imap #-}

instance FoldableWithIndex [Int] (Tree e) where
  ifoldMap f (Node _ lb ts) = f [] lb `mappend` ifoldMap (\i -> ifoldMap (f . (:) i)) ts
  {-# INLINE ifoldMap #-}

instance TraversableWithIndex [Int] (Tree e) where
  itraverse f (Node br lb ts) = Node br <$> f [] lb <*> itraverse (\i -> itraverse (f . (:) i)) ts
  {-# INLINE itraverse #-}

parseFileWith :: Parser a -> FilePath -> IO a
parseFileWith p f = do
  l <- if "gz" `isSuffixOf` f
    then BL.toStrict . decompress <$> BL.readFile f
    else BS.readFile f
  return $ either error id $ parseOnly p l

-- | Parse first Newick tree in file.
oneTree :: FilePath -> IO (Tree Double BS.ByteString)
oneTree f = do
  pt <- parseFileWith (newick Standard) f
  return $ first fromLength $ either error id $ phyloToLengthTree pt

-- | Parse one or more Newick trees until end of file.
someTrees :: FilePath -> IO [Tree Double BS.ByteString]
someTrees f = do
  pts <- parseFileWith (someNewick Standard) f
  return $ map (first fromLength . either error id . phyloToLengthTree) pts
