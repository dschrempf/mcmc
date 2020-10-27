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
    rootLabel,
    rootBranch,
  )
where

import Control.Lens
import ELynx.Tree

-- | Lens to a specific node.
subTreeAt :: Path -> Lens' (Tree e a) (Tree e a)
subTreeAt p =
  lens
    (getSubTreeUnsafe p)
    (\t t' -> let pos = goPathUnsafe p $ fromTree t in toTree $ pos {current = t'})

-- | Lens to the label of the root node.
rootLabel :: Lens' (Tree e a) a
rootLabel = lens label (\(Node br _ ts) lb -> Node br lb ts)

-- | Lens to the branch of the root node.
rootBranch :: Lens' (Tree e a) e
rootBranch = lens branch (\(Node _ lb ts) br -> Node br lb ts)
