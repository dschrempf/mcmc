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
  ( Path,
    root,
    mrca,
    constrainHard,
    constrainSoft,
    calibrate,
    calibrateUniform,
    calibrateUniformSoft,
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
isAncestor xs t = not $ any (`S.notMember` lvs) xs
  where
    lvs = S.fromList $ leaves t

isMrca :: Ord a => [a] -> Tree e a -> Bool
isMrca xs t = isAncestor xs t && all (not . isAncestor xs) (forest t)

-- | Path from the root of a tree to the node of the tree.
--
-- The position is specific to a tree topology. If the topology changes, the
-- position becomes invalid.
type Path = [Int]

-- | The position of the root.
root :: Path
root = []

-- | Get the path to the MRCA of the given list of node labels on the given
-- tree.
mrca :: Ord a => [a] -> Tree e a -> Maybe Path
mrca xs = (return . tail) <=< go 0
  where
    go i t
      | isMrca xs t = Just [i]
      | isAncestor xs t = Just $ i : head (catMaybes [go j t' | (j, t') <- zip [0 ..] (forest t)])
      | otherwise = Nothing

getHeightFromNode :: [Int] -> Tree Double Double -> Double
getHeightFromNode p = label . current . unsafeGoPath p . fromTree

-- | Hard constrain order of nodes with given paths using a truncated uniform
-- distribution.
--
-- Assume the branch and node labels denote branch length and node height,
-- respecitvely.
constrainHard ::
  -- | Young node (closer to the leaves).
  Path ->
  -- | Old node (closer to the root).
  Path ->
  Tree Double Double ->
  Log Double
constrainHard y o t
  | y `isPrefixOf` o = error "constrain: Young node is direct ancestor of old node (?)."
  | o `isPrefixOf` y = error "constrain: No need to constrain old node which is direct ancestor of young node."
  | getHeightFromNode y t < getHeightFromNode o t = 1
  | otherwise = 0

-- | Soft constrain order of nodes with given paths.
--
-- - When the node order is correct, a uniform distribution is used.
--
-- - When the node order is incorrect, a one-sided normal distribution with
--   given standard deviation is used. The normal distribution is normalized
--   such that the complete distribution of the constrained is continuous. Use
--   of the normal distribution also ensures that the first derivative is
--   continuous.
--
-- Assume the branch and node labels denote branch length and node height,
-- respecitvely.
constrainSoft ::
  -- | Standard deviation of one sided normal distribution.
  Double ->
  -- | Young node (closer to the leaves).
  Path ->
  -- | Old node (closer to the root).
  Path ->
  Tree Double Double ->
  Log Double
constrainSoft s y o t
  | y `isPrefixOf` o = error "constrain: Young node is direct ancestor of old node (?)."
  | o `isPrefixOf` y = error "constrain: No need to constrain old node which is direct ancestor of young node."
  | hY < hO = 1
  | otherwise = Exp $ logDensity d (hY - hO) - logDensity d 0
  where hY = getHeightFromNode y t
        hO = getHeightFromNode o t
        d = normalDistr 0 s

-- | Calibrate height of a node with given path using the normal distribution.
--
-- Assume the branch and node labels denote branch length and node height,
-- respecitvely.
calibrate ::
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  Path ->
  Tree Double Double ->
  Log Double
calibrate m s p = Exp . logDensity (normalDistr m s) . getHeightFromNode p

-- | Calibrate height of a node with given path using the uniform distribution.
--
-- Assume the branch and node labels denote branch length and node height,
-- respecitvely.
calibrateUniform ::
  -- | Lower bound.
  Double ->
  -- | Upper bound.
  Double ->
  Path ->
  Tree Double Double ->
  Log Double
calibrateUniform a b p t | h > b = 0
                         | h < a = 0
                         | otherwise = 1
  where h = getHeightFromNode p t

-- | Calibrate height of a node with given path.
--
-- - When the node is in the given bounds, a uniform distribution is used.
--
-- - When the node is out of bounds, a one-sided normal distribution with given
--   standard deviation is used. The normal distribution is normalized such that
--   the complete distribution of the constrained is continuous. Use of the
--   normal distribution also ensures that the first derivative is continuous.
--
-- Assume the branch and node labels denote branch length and node height,
-- respecitvely.
calibrateUniformSoft ::
  -- | Standard deviation of one sided normal distributions.
  Double ->
  -- | Lower bound.
  Double ->
  -- | Upper bound.
  Double ->
  Path ->
  Tree Double Double ->
  Log Double
calibrateUniformSoft s a b p t | h > b = Exp $ logDensity d (h - b) - logDensity d 0
                               | h < a = Exp $ logDensity d (a - h) - logDensity d 0
                               | otherwise = 1
  where h = getHeightFromNode p t
        d = normalDistr 0 s
