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
  ( NodePos,
    mrca,
    constrain,
    calibrate,
  )
where

import Data.List
import Data.Maybe
import Control.Monad
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

-- | The position of a node on a tree.
--
-- The position is specific to a tree topology. If the topology changes, the
-- position becomes invalid.
newtype NodePos = NP Path

-- | Get the path to the MRCA of the given list of node labels on the given
-- tree.
mrca :: Ord a => [a] -> Tree e a -> Maybe NodePos
mrca xs = (return . NP) <=< go 0
  where
    go i t
      | isMrca xs t = Just []
      | isAncestor xs t = Just $ i : head (catMaybes [go j t' | (j, t') <- zip [0 ..] (forest t)])
      | otherwise = Nothing

getHeight :: Measurable e => [Int] -> Tree e a -> Double
getHeight p = height . current . unsafeGoPath p . fromTree

-- | Constrain order of nodes with given paths using a truncated uniform
-- distribution.
constrain ::
  Measurable e =>
  -- | Young node (closer to the leaves).
  NodePos ->
  -- | Old node (closer to the root).
  NodePos ->
  Tree e a ->
  Log Double
constrain (NP y) (NP o) t
  | y `isPrefixOf` o = error "constrain: Young node is direct ancestor of old node (?)."
  | o `isPrefixOf` y = error "constrain: No need to constrain old node which is direct ancestor of young node."
  | getHeight y t < getHeight o t = 1
  | otherwise = 0

-- | Calibrate height of a node with given path using the normal distribution.
calibrate ::
  Measurable e =>
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  NodePos ->
  Tree e a ->
  Log Double
calibrate m s (NP p) = Exp . logDensity (normalDistr m s) . getHeight p
