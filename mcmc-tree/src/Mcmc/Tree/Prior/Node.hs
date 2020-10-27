-- |
-- Module      :  Mcmc.Tree.Prior.Node
-- Description :  Relative node order constraints and node calibrations
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jul 27 10:49:11 2020.
module Mcmc.Tree.Prior.Node
  ( constrainHard,
    constrainSoft,
    calibrate,
    calibrateUniform,
    calibrateUniformSoft,
  )
where

import Data.List
import ELynx.Tree
import Mcmc.Tree.Height
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal

-- Get the height of the node at path on the tree.
--
-- __Assume the node labels denote node height__.
getHeightFromNode :: HasHeight a => Path -> Tree Double a -> Double
getHeightFromNode p = getHeight . label . getSubTreeUnsafe p

-- | Hard constrain order of nodes with given paths using a truncated uniform
-- distribution.
--
-- __Assume the node labels denote node height__.
constrainHard ::
  HasHeight a =>
  -- | Path to younger node (closer to the leaves).
  Path ->
  -- | Path to older node (closer to the root).
  Path ->
  Tree Double a ->
  Log Double
constrainHard y o t
  | y `isPrefixOf` o = error "constrainHard: Young node is direct ancestor of old node (?)."
  | o `isPrefixOf` y = error "constrainHard: No need to constrain old node which is direct ancestor of young node."
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
-- __Assume the node labels denote node height__.
constrainSoft ::
  HasHeight a =>
  -- | Standard deviation of one sided normal distribution.
  Double ->
  -- | Path to younger node (closer to the leaves).
  Path ->
  -- | Path to older node (closer to the root).
  Path ->
  Tree Double a ->
  Log Double
constrainSoft s y o t
  | y `isPrefixOf` o = error "constrainSoft: Young node is direct ancestor of old node (?)."
  | o `isPrefixOf` y = error "constrainSoft: No need to constrain old node which is direct ancestor of young node."
  | hY < hO = 1
  | otherwise = Exp $ logDensity d (hY - hO) - logDensity d 0
  where
    hY = getHeightFromNode y t
    hO = getHeightFromNode o t
    d = normalDistr 0 s

-- | Calibrate height of a node with given path using the normal distribution.
--
-- __Assume the node labels denote node height__.
calibrate ::
  HasHeight a =>
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  Path ->
  Tree Double a ->
  Log Double
calibrate m s p = Exp . logDensity (normalDistr m s) . getHeightFromNode p

-- | Calibrate height of a node with given path using the uniform distribution.
--
-- __Assume the node labels denote node height__.
calibrateUniform ::
  HasHeight a =>
  -- | Lower bound.
  Double ->
  -- | Upper bound.
  Double ->
  Path ->
  Tree Double a ->
  Log Double
calibrateUniform a b p t
  | h > b = 0
  | h < a = 0
  | otherwise = 1
  where
    h = getHeightFromNode p t

-- | Calibrate height of a node with given path.
--
-- - When the node is in the given bounds, a uniform distribution is used.
--
-- - When the node is out of bounds, a one-sided normal distribution with given
--   standard deviation is used. The normal distribution is normalized such that
--   the complete distribution of the constrained is continuous. Use of the
--   normal distribution also ensures that the first derivative is continuous.
--
-- __Assume the node labels denote node height__.
calibrateUniformSoft ::
  HasHeight a =>
  -- | Standard deviation of one sided normal distributions.
  Double ->
  -- | Lower bound.
  Double ->
  -- | Upper bound.
  Double ->
  Path ->
  Tree Double a ->
  Log Double
calibrateUniformSoft s a b p t
  | h > b = Exp $ logDensity d (h - b) - logDensity d 0
  | h < a = Exp $ logDensity d (a - h) - logDensity d 0
  | otherwise = 1
  where
    h = getHeightFromNode p t
    d = normalDistr 0 s
