{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Tree.Lens
-- Description :  Lenses on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Aug 19 08:55:42 2020.
module Mcmc.Tree.Lens
  ( Path,
    subTreeAt,
    root,
    stem,
    lengthE,
  )
where

import Control.Lens
import ELynx.Tree

splitAt' :: Int -> [a] -> ([a], a, [a])
splitAt' i xs = (ls, head rs, tail rs)
  where
    (ls, rs) = splitAt i xs

assemble :: e -> a -> [Tree e a] -> [Tree e a] -> Tree e a -> Tree e a
assemble br lb ls rs c = Node br lb $ ls ++ (c : rs)

-- | A specific sub tree.
subTreeAt :: Path -> Lens' (Tree e a) (Tree e a)
subTreeAt pth f s = go s pth
  where
    go t [] = f t
    go (Node lb br ts) (x : xs) =
      let (ls, c, rs) = splitAt' x ts
       in assemble lb br ls rs <$> go c xs

-- -- Around 10 percent slower for trees with five to ten levels, because they
-- -- have to be traversed twice. However, the loss of speed may be worse for
-- -- deeper structures.
--
-- -- | A specific sub tree.
-- subTreeAt :: Path -> Lens' (Tree e a) (Tree e a)
-- subTreeAt p =
--   lens
--     (getSubTreeUnsafe p)
--     (\t t' -> let pos = goPathUnsafe p $ fromTree t in toTree $ pos {current = t'})

-- | Label of the root node.
root :: Lens' (Tree e a) a
root = lens label (\(Node br _ ts) lb -> Node br lb ts)

-- | Branch attached to the root node.
stem :: Lens' (Tree e a) e
stem = lens branch (\(Node _ lb ts) br -> Node br lb ts)

-- | Length. Setter calls 'error' if length is negative.
lengthE:: Lens' Length Double
lengthE = lens fromLength (\_ x -> either error id $ toLength x)
