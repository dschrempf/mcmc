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
  ( -- * Constraints
    constrainHard,
    constrainSoft,

    -- * Calibtrations
    NonNegative,
    ExtendedPositive,
    Interval,
    properInterval,
    lowerBoundOnly,
    transformInterval,
    calibrate,
    calibrateUniform,
    calibrateUniformSoft,
  )
where

import Control.Lens
import Data.List
import ELynx.Tree
import Mcmc.Tree.Lens
import Mcmc.Tree.Types
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import Text.Read

-- Get the height of the node at path on the tree.
--
-- __Assume the node labels denote node height__.
getHeightFromNode :: HasHeight a => Path -> Tree e a -> Height
getHeightFromNode p t = t ^. subTreeAtUnsafeL p . labelL . hasHeightL

-- | Hard constrain order of nodes with given paths using a truncated uniform
-- distribution.
constrainHard ::
  HasHeight a =>
  -- | Path to younger node (closer to the leaves).
  Path ->
  -- | Path to older node (closer to the root).
  Path ->
  Tree e a ->
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
constrainSoft ::
  HasHeight a =>
  -- | Standard deviation of one sided normal distribution.
  Double ->
  -- | Path to younger node (closer to the leaves).
  Path ->
  -- | Path to older node (closer to the root).
  Path ->
  Tree e a ->
  Log Double
constrainSoft s y o t
  | y `isPrefixOf` o = error "constrainSoft: Young node is direct ancestor of old node (?)."
  | o `isPrefixOf` y = error "constrainSoft: No need to constrain old node which is direct ancestor of young node."
  | hY < hO = 1
  | otherwise = Exp $ logDensity d (hY - hO) - logDensity d 0
  where
    hY = fromHeight $ getHeightFromNode y t
    hO = fromHeight $ getHeightFromNode o t
    d = normalDistr 0 s

-- | Non-negative number.
newtype NonNegative = NonNegative {fromNonNegative :: Double}
  deriving (Eq)

nonNegative :: Double -> NonNegative
nonNegative x
  | x < 0 = error "nonNegative: Negative value."
  | otherwise = NonNegative x

instance Read NonNegative where
  readPrec = do
    x <- readPrec
    return $ nonNegative x

instance Show NonNegative where
  showsPrec p (NonNegative x) = showsPrec p x

-- | Positive number or infinity.
data ExtendedPositive = Positive Double | Infinity
  deriving (Eq)

positive :: Double -> ExtendedPositive
positive x
  | x <= 0 = error "positive: Zero or negative value."
  | otherwise = Positive x

positiveReadPrec :: ReadPrec ExtendedPositive
positiveReadPrec = do
  x <- readPrec
  return $ positive x

infinityReadPrec :: ReadPrec ExtendedPositive
infinityReadPrec = do
  Ident "Infinity" <- lexP
  return Infinity

instance Read ExtendedPositive where
  readPrec = positiveReadPrec <++ infinityReadPrec

instance Show ExtendedPositive where
  showsPrec p (Positive x) = showsPrec p x
  showsPrec p Infinity = showsPrec p "Infinity"

-- | Open interval \((a,b)\) with \(a < b\), \(a \in [0, \infty)\) and \(b \in
-- (0, \infty]\).
data Interval = Interval NonNegative ExtendedPositive
  deriving (Eq)

instance Show Interval where
  show (Interval a b) = "(" ++ show a ++ ", " ++ show b ++ ")"

-- | Specify a lower and an upper bound.
properInterval :: Double -> Double -> Interval
properInterval a b
  | a < b = Interval (nonNegative a) (positive b)
  | otherwise = error "properInterval: Left bound equal or larger right bound."

-- | Specify a lower bound only. The upper bound is set to infinity.
lowerBoundOnly :: Double -> Interval
lowerBoundOnly a = Interval (nonNegative a) Infinity

-- | Transform an interval by applying a multiplicative change.
--
-- Useful when the tree is normalized and height values have to be converted
-- from relative heights to absolute heights.
transformInterval :: Double -> Interval -> Interval
transformInterval x (Interval a b)
  | x <= 0 = error "transform: Multiplier is zero or negative."
  | otherwise = Interval a' b'
  where
    a' = NonNegative $ x * fromNonNegative a
    b' = case b of
      Positive y -> Positive $ x * y
      Infinity -> Infinity

-- No number is bigger than a non-existing upper bound..
(>*) :: Double -> ExtendedPositive -> Bool
_ >* Infinity = False
h >* Positive b = h > b

-- | Calibrate height of a node with given path using the normal distribution.
calibrate ::
  HasHeight a =>
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  Path ->
  Tree e a ->
  Log Double
calibrate m s p = Exp . logDensity (normalDistr m s) . fromHeight . getHeightFromNode p

-- | Calibrate height of a node with given path using the uniform distribution.
--
-- If the upper bound is not given, no upper bound is used.
calibrateUniform ::
  HasHeight a =>
  Interval ->
  Path ->
  Tree e a ->
  Log Double
calibrateUniform (Interval a b) p t
  | h <= a' = 0
  | h >* b = 0
  | otherwise = 1
  where
    a' = fromNonNegative a
    h = fromHeight $ getHeightFromNode p t

-- | Calibrate height of a node with given path.
--
-- - When the node is in the given bounds, a uniform distribution is used.
--
-- - When the node is out of bounds, a one-sided normal distribution with given
--   standard deviation is used. The normal distribution is normalized such that
--   the complete distribution of the constrained is continuous. Use of the
--   normal distribution also ensures that the first derivative is continuous.
--
-- If the upper bound is not given, no upper bound is used.
calibrateUniformSoft ::
  HasHeight a =>
  -- | Standard deviation of one sided normal distributions.
  Double ->
  Interval ->
  Path ->
  Tree e a ->
  Log Double
calibrateUniformSoft s (Interval a b) p t
  | h <= a' = Exp $ logDensity d (a' - h) - logDensity d 0
  | h >* b = case b of
    Infinity -> 1.0
    Positive b' -> Exp $ logDensity d (h - b') - logDensity d 0
  | otherwise = 1
  where
    a' = fromNonNegative a
    h = fromHeight $ getHeightFromNode p t
    d = normalDistr 0 s
