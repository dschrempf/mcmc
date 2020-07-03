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
-- Creation date: Fri Jul  3 09:35:43 2020.
--
-- Some functions are inspired by
-- [Biobase.Newick.Import](https://hackage.haskell.org/package/BiobaseNewick).
--
-- [Specifications](http://evolution.genetics.washington.edu/phylip/newicktree.html)
--
-- - In particular, no conversion from _ to (space) is done right now.
module Mcmc.Tree.Newick
  ( Parser,
    newick,
    oneNewick,
    manyNewick,
  )
where

import Algebra.Graph.Labelled.AdjacencyMap
import qualified Data.ByteString.Lazy as L
import Data.ByteString.Lazy (ByteString)
import Data.ByteString.Internal (c2w)
import Data.Void
import Data.Word
import Mcmc.Tree.Tree
import Text.Megaparsec
import Text.Megaparsec.Byte
import Text.Megaparsec.Byte.Lexer
  ( decimal,
    float,
  )

-- | Parse a single Newick tree. Also succeeds when more trees follow.
newick :: ByteString -> Tree Double ByteString
newick = parseByteStringWith "newick" newickP

-- | Parse a single Newick tree. Fails when end of file is not reached.
oneNewick :: ByteString -> Tree Double ByteString
oneNewick = parseByteStringWith "oneNewick" oneNewickP

-- | Parse many Newick trees until end of file.
manyNewick :: ByteString -> [Tree Double ByteString]
manyNewick = parseByteStringWith "manyNewick" manyNewickP

-- Parse a 'ByteString' and extract the result.
parseByteStringWith
  :: (ShowErrorComponent e)
  => String                  -- ^ Name of byte string.
  -> Parsec e ByteString a -- ^ Parser.
  -> ByteString            -- ^ Input.
  -> a
parseByteStringWith s p l = case parse p s l of
  Left  err -> error $ errorBundlePretty err
  Right val -> val

type Parser = Parsec Void ByteString

w :: Char -> Parser Word8
w = char . c2w

-- Parse a single Newick tree. Also succeeds when more trees follow.
newickP :: Parser (Tree Double ByteString)
newickP = space *> tree <* w ';' <* space <?> "newick"

-- Parse a single Newick tree. Fails when end of file is not reached.
oneNewickP :: Parser (Tree Double ByteString)
oneNewickP = newickP <* eof <?> "oneNewick"

-- Parse many Newick trees until end of file.
manyNewickP :: Parser [Tree Double ByteString]
manyNewickP = some newickP <* eof <?> "manyNewick"

tree :: Parser (Tree Double ByteString)
tree = branched <|> leaf <?> "tree"

branched :: Parser (Tree Double ByteString)
branched = do
  f <- fmap fromTree <$> forest
  n <- node
  b <- toD <$> branchLength <?> "branched"
  return $ Tree $ connect b (vertex n) (overlays f)

-- A forest is a set of trees separated by @,@ and enclosed by parentheses.
forest :: Parser [Tree Double ByteString]
forest = between (w '(') (w ')') (tree `sepBy1` w ',') <?> "forest"

-- A leaf is a node without children.
leaf :: Parser (Tree Double ByteString)
leaf = do
  n <- node <?> "leaf"
  return $ Tree $ vertex n

-- A node name can be any string of printable characters except blanks, colons,
-- semicolons, parentheses, and square brackets (and commas).
node :: Parser ByteString
node = L.pack <$> many (satisfy checkNameCharacter) <?> "node"

checkNameCharacter :: Word8 -> Bool
checkNameCharacter c = c `notElem` map c2w " :;()[],"

-- Branch length.
branchLength :: Parser Double
branchLength = w ':' *> branchLengthGiven <?> "branchLength"

branchLengthGiven :: Parser Double
branchLengthGiven = try float <|> decimalAsDouble

decimalAsDouble :: Parser Double
decimalAsDouble = fromIntegral <$> (decimal :: Parser Int)
