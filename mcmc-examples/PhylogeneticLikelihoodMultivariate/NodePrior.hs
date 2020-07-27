-- |
-- Module      :  NodePrior
-- Description :  Relative node order constraints and node calibrations
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jul 27 10:49:11 2020.
module NodePrior
  ( getPathToMrca,
    constrain,
    calibrate,
  )
where

import Data.Maybe
import qualified Data.Set as S
import ELynx.Data.Tree
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal

isAncestor :: Ord a => [a] -> Tree e a -> Bool
isAncestor xs t = any (`S.notMember` lvs) xs
  where
    lvs = S.fromList $ leaves t

isMrca :: Ord a => [a] -> Tree e a -> Bool
isMrca xs t = isAncestor xs t && all (not . isAncestor xs) (forest t)

type Path = [Int]

getPathToMrca :: Ord a => [a] -> Tree e a -> Maybe [Int]
getPathToMrca = go 0
  where
    go i xs t
      | isMrca xs t = Just []
      | isAncestor xs t = Just $ i : head (catMaybes [go j xs t' | (j, t') <- zip [0 ..] (forest t)])
      | otherwise = Nothing

getHeight :: Measurable e => Path -> Tree e a -> Double
getHeight p = height . current . unsafeGoPath p . fromTree

constrain :: Measurable e => Path -> Path -> Tree e a -> Log Double
constrain p1 p2 t
  | getHeight p1 t < getHeight p2 t = 1
  | otherwise = 0

-- TODO: This only honors the tree branch lengths, but not a global modifier.
calibrate :: Measurable e => Path -> Tree e a -> Log Double
calibrate p = Exp . logDensity standard . getHeight p
