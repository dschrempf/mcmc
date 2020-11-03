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
    lengthError,
    lengthUnsafe,
  )
where

import Control.Lens
import ELynx.Tree

splitAt' :: Int -> [a] -> ([a], a, [a])
splitAt' i xs = (ls, head rs, tail rs)
  where
    (ls, rs) = splitAt i xs
{-# INLINE splitAt' #-}

-- type Lens s t a b = forall f. Functor f => (a -> f b) -> s -> f t
-- type Lens' s a = forall f. Functor f => (a -> f a) -> s -> f s
-- lens :: (s -> a) -> (s -> b -> t) -> Lens s t a b
-- lens sa sbt afb s = sbt s <$> afb (sa s)

-- | A specific sub tree.
--
-- Call 'error' if the path is invalid.
subTreeAt :: Path -> Lens' (Tree e a) (Tree e a)
subTreeAt pth f s = go s pth
  where
    go t [] = f t
    go (Node lb br ts) (x : xs) =
      let (ls, c, rs) = splitAt' x ts
       in assemble lb br ls rs <$> go c xs
    assemble :: e -> a -> Forest e a -> Forest e a -> Tree e a -> Tree e a
    assemble br lb ls rs c = Node br lb $ ls ++ (c : rs)

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
root f (Node br lb ts) = assemble br ts <$> f lb
  where
    assemble :: e -> Forest e a -> a -> Tree e a
    assemble br' ts' lb' = Node br' lb' ts'

-- | Branch attached to the root node.
stem :: Lens' (Tree e a) e
stem f (Node br lb ts) = assemble lb ts <$> f br
  where
    assemble :: a -> Forest e a -> e -> Tree e a
    assemble lb' ts' br' = Node br' lb' ts'

-- | Length.
--
-- Call 'error' if length is negative.
lengthError :: Lens' Length Double
lengthError f l = either error id . toLength <$> f (fromLength l)


-- | Length.
--
-- Non-negativity property is not ensured.
lengthUnsafe :: Lens' Length Double
lengthUnsafe f l = toLengthUnsafe <$> f (fromLength l)
