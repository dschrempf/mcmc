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
    Constraint (..),
    constraint,
    validateConstraint,
    constrainHardUnsafe,
    constrainSoftUnsafe,

    -- * Intervals
    NonNegative,
    ExtendedPositive,
    Interval,
    properInterval,
    lowerBoundOnly,
    transformInterval,

    -- * Calibrations
    Calibration (..),
    calibration,
    calibrateHardUnsafe,
    calibrateSoftUnsafe,
  )
where

import Control.Lens
import Data.List
import ELynx.Tree
import Mcmc.Tree.Lens
import Mcmc.Tree.Mrca
import Mcmc.Tree.Types
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import Text.Read

-- | Constraints define node orders.
--
-- For example,
--
-- @
--   Constraint "Name" YOUNGER OLDER
-- @
--
-- ensures that the node with path @YOUNGER@ is younger than the node with path
-- @OLDER@.
data Constraint = Constraint
  { constraintName :: String,
    -- | Path to younger node (closer to the leaves).
    constraintYoungNode :: Path,
    -- | Path to older node (closer to the root).
    constraintOldNode :: Path
  }
  deriving (Eq, Read, Show)

-- Get the height of the node at path on the tree.
--
-- __Assume the node labels denote node height__.
getHeightFromNodeUnsafe :: HasHeight a => Path -> Tree e a -> Height
getHeightFromNodeUnsafe p t = t ^. subTreeAtUnsafeL p . labelL . hasHeightL

-- | Check if a constraint is valid.
--
-- Return 'Left' if:
--
-- - The younger node is a direct ancestor of the old node.
--
-- - The older node is a direct ancestor of the young node.
validateConstraint ::
  Constraint ->
  Either String Constraint
validateConstraint c
  | y `isPrefixOf` o = Left $ getErrMsg "Younger node is direct ancestor of older node (?)."
  | o `isPrefixOf` y = Left $ getErrMsg "No need to constrain older node which is direct ancestor of younger node."
  | otherwise = Right c
  where
    y = constraintYoungNode c
    o = constraintOldNode c
    getErrMsg msg = "validateConstraint: " ++ constraintName c ++ ": " ++ msg

-- | Create and validate a constraint.
--
-- Call 'error' if:
--
-- - A node cannot be found on the tree.
--
-- - The younger node is a direct ancestor of the old node.
--
-- - The older node is a direct ancestor of the young node.
constraint ::
  (Ord a, Show a) =>
  Tree e a ->
  -- | Name.
  String ->
  -- | The most recent common ancestor of the given leaves is the younger node.
  [a] ->
  -- | The most recent common ancestor of the given leave is the older node.
  [a] ->
  Constraint
constraint t n ys os =
  either error id $
    validateConstraint $
      Constraint n y o
  where
    err msg = error $ "constraint: " ++ n ++ ": " ++ msg
    y = either err id $ mrca ys t
    o = either err id $ mrca os t

-- | Hard constrain order of nodes with given paths.
--
-- A truncated, improper uniform distribution is used.
--
-- For reasons of computational efficiency, the paths are not checked for
-- validity. Please do so beforehand using 'constraint' or 'validateConstraint'.
constrainHardUnsafe ::
  HasHeight a =>
  Constraint ->
  Tree e a ->
  Log Double
constrainHardUnsafe c t
  | getHeightFromNodeUnsafe y t < getHeightFromNodeUnsafe o t = 1
  | otherwise = 0
  where
    y = constraintYoungNode c
    o = constraintOldNode c

-- | Soft constrain order of nodes with given paths.
--
-- When the node order is correct, a uniform distribution is used.
--
-- When the node order is incorrect, a one-sided normal distribution with given
-- standard deviation is used. The normal distribution is normalized such that
-- the complete distribution of the constraint is continuous. Use of the normal
-- distribution also ensures that the first derivative is continuous.
--
-- For reasons of computational efficiency, the paths are not checked for
-- validity. Please do so beforehand using 'constraint' or 'validateConstraint'.
constrainSoftUnsafe ::
  HasHeight a =>
  -- | Standard deviation of one sided normal distribution.
  Double ->
  Constraint ->
  Tree e a ->
  Log Double
constrainSoftUnsafe s c t
  | hY < hO = 1
  | otherwise = Exp $ logDensity d (hY - hO) - logDensity d 0
  where
    hY = fromHeight $ getHeightFromNodeUnsafe y t
    hO = fromHeight $ getHeightFromNodeUnsafe o t
    d = normalDistr 0 s
    y = constraintYoungNode c
    o = constraintOldNode c

-- | Non-negative number.
newtype NonNegative = NonNegative {fromNonNegative :: Double}
  deriving (Eq)

nonNegative :: Double -> NonNegative
nonNegative x
  | x < 0 = error "nonNegative: Negative value."
  | otherwise = NonNegative x

instance Read NonNegative where
  readPrec = nonNegative <$> readPrec

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
positiveReadPrec = positive <$> readPrec

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

-- | A calibration is specified by a name, a node at given path, and height
-- boundaries.
--
-- For example,
--
-- @
--   let c = Calibration "Root" [] YOUNG OLD
-- @
--
-- ensures that the root node is older than @YOUNG@, and younger than @OLD@.
data Calibration = Calibration
  { calibrationName :: String,
    calibrationNode :: Path,
    calibrationInterval :: Interval
  }
  deriving (Eq, Show)

-- | Create and validate a calibration.
--
-- Call 'error' if:
--
-- - The node cannot be found on the tree.
calibration ::
  (Ord a, Show a) =>
  Tree e a ->
  -- | Name.
  String ->
  -- | The most recent common ancestor of the given leaves is the calibrated node.
  [a] ->
  Interval ->
  Calibration
calibration t n xs = Calibration n p
  where
    err msg = error $ "calibration: " ++ n ++ ": " ++ msg
    p = either err id $ mrca xs t

-- | Calibrate height of a node with given path using the uniform distribution.
--
-- If the upper bound is not given, no upper bound is used.
--
-- For reasons of computational efficiency, the path is not checked for
-- validity. Please do so beforehand using 'calibration'.
calibrateHardUnsafe ::
  HasHeight a =>
  Calibration ->
  Tree e a ->
  Log Double
calibrateHardUnsafe c t
  | h <= a' = 0
  | h >* b = 0
  | otherwise = 1
  where
    a' = fromNonNegative a
    h = fromHeight $ getHeightFromNodeUnsafe p t
    (Interval a b) = calibrationInterval c
    p = calibrationNode c

-- | Calibrate height of a node with given path.
--
-- When the node is in the given bounds, a uniform distribution is used.
--
-- When the node is out of bounds, a one-sided normal distribution with given
-- standard deviation is used. The normal distribution is normalized such that
-- the complete distribution of the constrained is continuous. Use of the normal
-- distribution also ensures that the first derivative is continuous.
--
-- If the upper bound is not given, no upper bound is used.
--
-- For reasons of computational efficiency, the path is not checked for
-- validity. Please do so beforehand using 'calibration'.
calibrateSoftUnsafe ::
  HasHeight a =>
  -- | Standard deviation of one sided normal distributions.
  Double ->
  Calibration ->
  Tree e a ->
  Log Double
calibrateSoftUnsafe s c t
  | h <= a' = Exp $ logDensity d (a' - h) - logDensity d 0
  | h >* b = case b of
    Infinity -> 1.0
    Positive b' -> Exp $ logDensity d (h - b') - logDensity d 0
  | otherwise = 1
  where
    a' = fromNonNegative a
    h = fromHeight $ getHeightFromNodeUnsafe p t
    d = normalDistr 0 s
    (Interval a b) = calibrationInterval c
    p = calibrationNode c
