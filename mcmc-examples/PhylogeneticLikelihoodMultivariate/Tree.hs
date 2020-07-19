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
  ( oneTree,
    someTrees,
  )
where

import Codec.Compression.GZip
import Data.ByteString.Lazy (ByteString)
import qualified Data.ByteString.Lazy as L
import ELynx.Data.Tree
import ELynx.Import.Tree.Newick
import Text.Megaparsec

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
someTrees :: FilePath -> IO [Tree Length Int]
someTrees p = do
  pts <- parseFileWith "someNewick" (someNewick Standard) p
  return $ map (either error id . phyloToLengthTree . identify) pts
