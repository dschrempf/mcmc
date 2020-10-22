{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Tools
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Oct 22 16:53:59 2020.
module Tools
  ( getBranches,
    sumFirstTwo,
  )
where

import qualified Data.Vector.Storable as V
import ELynx.Tree

-- | Get all branches of a rooted tree. Store the branches in a vector such that
-- the two branches leading to the root are the first two entries of the vector.
-- Ignore the root branch.
getBranches :: Tree Double a -> V.Vector Double
getBranches (Node _ _ [l, r]) = V.fromList $ head ls : head rs : tail ls ++ tail rs
  where
    ls = branches l
    rs = branches r
getBranches _ = error "getBranches: Root node is not bifurcating."

-- | Sum the first two elements of a vector. Needed to merge the two branches
-- leading to the root.
sumFirstTwo :: V.Vector Double -> V.Vector Double
sumFirstTwo v = (v V.! 0 + v V.! 1) `V.cons` V.drop 2 v
