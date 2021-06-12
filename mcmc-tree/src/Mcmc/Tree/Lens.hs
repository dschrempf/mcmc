{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Tree.Lens
-- Description :  Lenses on trees
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Aug 19 08:55:42 2020.
module Mcmc.Tree.Lens
  ( Path,

    -- ** Tree lenses
    labelL,
    forestL,
    branchL,
    subTreeAtUnsafeL,

    -- ** Tree zippers
    currentL,

    -- ** Height and length lenses
    hasHeightL,
    heightUnsafeL,
    hasLengthL,
    lengthUnsafeL,
  )
where

import Control.Lens
import ELynx.Tree
import Mcmc.Tree.Types

-- | Branch attached to the root node.
branchL :: Lens' (Tree e a) e
branchL f (Node br lb ts) = assemble lb ts <$> f br
  where
    assemble :: a -> Forest e a -> e -> Tree e a
    assemble lb' ts' br' = Node br' lb' ts'

-- | Label of the root node.
labelL :: Lens' (Tree e a) a
labelL f (Node br lb ts) = assemble br ts <$> f lb
  where
    assemble :: e -> Forest e a -> a -> Tree e a
    assemble br' ts' lb' = Node br' lb' ts'

-- | Forest of the root node.
forestL :: Lens' (Tree e a) (Forest e a)
forestL f (Node br lb ts) = Node br lb <$> f ts

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
subTreeAtUnsafeL :: Path -> Lens' (Tree e a) (Tree e a)
subTreeAtUnsafeL pth f s = go s pth
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
-- subTreeAtUnsafeL :: Path -> Lens' (Tree e a) (Tree e a)
-- subTreeAtUnsafeL p =
--   lens
--     (getSubTreeUnsafe p)
--     (\t t' -> let pos = goPathUnsafe p $ fromTree t in toTree $ pos {current = t'})

-- | Current focus of a zipper.
currentL :: Lens' (TreePos e a) (Tree e a)
currentL = lens current (\x t -> x {current = t})

-- | Height.
hasHeightL :: HasHeight a => Lens' a Height
hasHeightL = lens getHeight (flip setHeight)

-- | Height, unsafe.
--
-- Non-negativity is not ensured.
heightUnsafeL :: Lens' Height Double
heightUnsafeL f l = toHeightUnsafe <$> f (fromHeight l)

-- | Length of measurable types.
hasLengthL :: HasLength a => Lens' a Length
hasLengthL = lens getLen (flip setLen)

-- | Length, unsafe.
--
-- Non-negativity is not ensured.
lengthUnsafeL :: Lens' Length Double
lengthUnsafeL f l = toLengthUnsafe <$> f (fromLength l)
