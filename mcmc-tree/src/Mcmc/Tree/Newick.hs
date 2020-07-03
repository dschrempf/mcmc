-- |
-- Module      :  Mcmc.Tree.Newick
-- Description :  Convert to and from Newick format
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
--
-- Creation date: Fri Jul  3 09:35:43 2020.
--
-- Some functions are inspired by
-- [Biobase.Newick.Import](https://hackage.haskell.org/package/BiobaseNewick).
--
-- [Specifications](http://evolution.genetics.washington.edu/phylip/newicktree.html)
--
-- In particular, no conversion from _ to (space) is done right now.
module Mcmc.Tree.Newick
  ( newick,
    oneNewick,
    manyNewick,
  )
where

import Data.ByteString.Internal (c2w)
import qualified Data.ByteString.Lazy as L
import Data.ByteString.Lazy (ByteString)
import Data.Maybe
import Data.Tree
import Data.Void
import Data.Word
import Text.Megaparsec
import Text.Megaparsec.Byte
import Text.Megaparsec.Byte.Lexer
  ( decimal,
    float,
  )

-- | Parse a single Newick tree. Also succeeds when more trees follow.
newick :: FilePath -> IO (Tree (Double, ByteString))
newick = parseFileWith "newick" newickTree

-- | Parse a single Newick tree. Fails when end of file is not reached.
oneNewick :: FilePath -> IO (Tree (Double, ByteString))
oneNewick = parseFileWith "oneNewick" oneNewickTree

-- | Parse many Newick trees until end of file.
manyNewick :: FilePath -> IO [Tree (Double, ByteString)]
manyNewick = parseFileWith "manyNewick" manyNewickTree

parseFileWith ::
  (ShowErrorComponent e) =>
  -- | Name of byte string.
  String ->
  -- | Parser.
  Parsec e ByteString a ->
  -- | Input file.
  FilePath ->
  IO a
parseFileWith s p f = do
  l <- L.readFile f
  return $ case parse p s l of
    Left err -> error $ errorBundlePretty err
    Right val -> val

type L = (Double, ByteString)

type P = Parsec Void ByteString

w :: Char -> P Word8
w = char . c2w

-- Parse a single Newick tree. Also succeeds when more trees follow.
newickTree :: P (Tree L)
newickTree = space *> tree <* w ';' <* space <?> "newick"

-- Parse a single Newick tree. Fails when end of file is not reached.
oneNewickTree :: P (Tree L)
oneNewickTree = newickTree <* eof <?> "oneNewick"

-- Parse many Newick trees until end of file.
manyNewickTree :: P [Tree L]
manyNewickTree = some newickTree <* eof <?> "manyNewick"

tree :: P (Tree L)
tree = branched <|> leaf <?> "tree"

branched :: P (Tree L)
branched = do
  f <- forest
  l <- nodeLabel <?> "branched"
  return $ Node l f

-- A forest is a set of trees separated by @,@ and enclosed by parentheses.
forest :: P [Tree L]
forest = between (w '(') (w ')') (tree `sepBy1` w ',') <?> "forest"

-- A leaf is a node without children.
leaf :: P (Tree L)
leaf = do
  l <- nodeLabel <?> "leaf"
  return $ Node l []

-- A node has a name and a branch length.
nodeLabel :: P L
nodeLabel = do
  n <- name
  b <- branchLength
  return (fromMaybe 0 b, n)

checkNameCharacter :: Word8 -> Bool
checkNameCharacter c = c `notElem` map c2w " :;()[],"

-- A name can be any string of printable characters except blanks, colons,
-- semicolons, parentheses, and square brackets (and commas).
name :: P ByteString
name = L.pack <$> many (satisfy checkNameCharacter) <?> "name"

-- Branch length.
branchLength :: P (Maybe Double)
branchLength = optional $ w ':' *> (try float <|> decimalAsDouble)

decimalAsDouble :: P Double
decimalAsDouble = fromIntegral <$> (decimal :: P Int)
