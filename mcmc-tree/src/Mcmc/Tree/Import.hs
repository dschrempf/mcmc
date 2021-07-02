{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Mcmc.Tree.Import
-- Description :  Markov chain Monte Carlo sampling on trees
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri Jul  3 09:05:09 2020.
module Mcmc.Tree.Import
  ( oneTree,
    someTrees,
  )
where

import Codec.Compression.GZip
import Control.Lens
import Data.Attoparsec.Lazy
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List
import ELynx.Tree

instance FunctorWithIndex [Int] (Tree e) where
  imap f (Node br lb ts) = Node br (f [] lb) $ imap (\i -> imap (f . (:) i)) ts
  {-# INLINE imap #-}

instance FoldableWithIndex [Int] (Tree e) where
  ifoldMap f (Node _ lb ts) = f [] lb <> ifoldMap (\i -> ifoldMap (f . (:) i)) ts
  {-# INLINE ifoldMap #-}

instance TraversableWithIndex [Int] (Tree e) where
  itraverse f (Node br lb ts) = Node br <$> f [] lb <*> itraverse (\i -> itraverse (f . (:) i)) ts
  {-# INLINE itraverse #-}

parseFileWith :: Parser a -> FilePath -> IO a
parseFileWith p f = do
  l <-
    if "gz" `isSuffixOf` f
      then BL.toStrict . decompress <$> BL.readFile f
      else BS.readFile f
  return $ either error id $ parseOnly p l

-- | Parse first Newick tree in file.
oneTree :: NewickFormat -> FilePath -> IO (Tree Length Name)
oneTree fm f = do
  t <- parseFileWith (newick fm) f
  return $ either error id $ toLengthTree t

-- | Parse one or more Newick trees until end of file.
someTrees :: NewickFormat -> FilePath -> IO [Tree Length Name]
someTrees fm fn = do
  pts <- parseFileWith (someNewick fm) fn
  return $ map (either error id . toLengthTree) pts
