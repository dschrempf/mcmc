{-# LANGUAGE DeriveGeneric #-}

-- |
-- Module      :  Mcmc.Tree.Prior.Node.Calibration
-- Description :  Relative node order constraints and node calibrations
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jul 27 10:49:11 2020.
module Mcmc.Tree.Prior.Node.Calibration
  ( -- * Intervals
    NonNegative,
    ExtendedPositive,
    Interval,
    properInterval,
    lowerBoundOnly,
    transformInterval,

    -- * Calibrations
    Calibration (..),
    calibration,
    loadCalibrations,
    calibrateHard,
    calibrateSoft,
    calibrate,
  )
where

import Control.Monad
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Csv hiding (Name)
import Data.Function
import Data.List
import qualified Data.Vector as V
import ELynx.Tree hiding (partition)
import GHC.Generics
import Mcmc.Chain.Chain
import Mcmc.Statistics.Types
import Mcmc.Tree.Mrca
import Mcmc.Tree.Prior.Node.Common
import Mcmc.Tree.Types
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import Text.Read

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
  showsPrec _ Infinity = showString "Infinity"

-- | Open interval \((a,b)\) with \(a < b\), \(a \in [0, \infty)\) and \(b \in
-- (0, \infty]\).
data Interval = Interval NonNegative ExtendedPositive
  deriving (Eq)

instance Show Interval where
  show (Interval a b) = "(" ++ show a ++ ", " ++ show b ++ ")"

-- | Specify a lower and an upper bound.
properInterval :: LowerBoundary -> UpperBoundary -> Interval
properInterval a b
  | a < b = Interval (nonNegative a) (positive b)
  | otherwise = error "properInterval: Left bound equal or larger right bound."

-- | Specify a lower bound only. The upper bound is set to infinity.
lowerBoundOnly :: LowerBoundary -> Interval
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
      Positive bVal -> Positive $ x * bVal
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

prettyPrintCalibration :: Calibration -> String
prettyPrintCalibration (Calibration n p i) =
  "Calibration: "
    <> n
    <> " with path "
    <> show p
    <> " and interval "
    <> show i
    <> "."

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

data CalibrationData = CalibrationData String String String Double (Maybe Double)
  deriving (Generic, Show)

instance FromRecord CalibrationData

calibrationDataToCalibration :: Tree e Name -> CalibrationData -> Calibration
calibrationDataToCalibration t (CalibrationData n a b l mr) = calibration t n [a', b'] i
  where
    a' = Name $ BL.pack a
    b' = Name $ BL.pack b
    i = case mr of
      Nothing -> lowerBoundOnly l
      Just r -> properInterval l r

-- Get duplicate pairs of a list.
findDupsBy :: (a -> a -> Bool) -> [a] -> [[a]]
findDupsBy _ [] = []
findDupsBy eq (x : xs) = case partition (eq x) xs of
  ([], _) -> findDupsBy eq xs
  (ys, xs') -> (x : ys) : findDupsBy eq xs'

-- | Load and validate calibrations from file.
--
-- The calibration file is a comma separated values (CSV) file with rows of the
-- following format:
--
-- @
-- CalibrationName,LeafA,LeafB,LowerBoundary,UpperBoundary
-- @
--
-- The calibrated node is uniquely defined as the most recent common ancestor
-- (MRCA) of @LeafA@ and @LeafB@. The UpperBoundary can be omitted.
--
-- The following line defines a calibration with a lower boundary only:
--
-- @
-- Primates,Human,Chimpanzees,1e6,
-- @
--
-- Call 'error' if
--
-- - The file contains errors.
--
-- - An MRCA cannot be found.
--
-- - Redundant or conflicting calibrations are found (i.e., multiple
--   calibrations affect single nodes).
loadCalibrations :: Tree e Name -> FilePath -> IO (V.Vector Calibration)
loadCalibrations t f = do
  d <- BL.readFile f
  let mr = decode NoHeader d :: Either String (V.Vector CalibrationData)
      cds = either error id mr
  when (V.null cds) $ error $ "loadCalibrations: No calibrations found in file: " <> f <> "."
  let calsAll = V.map (calibrationDataToCalibration t) cds
  -- Check for redundant or conflicting calibrations.
  --
  let calsDupl = findDupsBy ((==) `on` calibrationNode) $ V.toList calsAll
  unless (null calsDupl) $ do
    -- Calibrations could also be removed. But then, which one should be removed?
    let render xs = unlines $ "Redundant and/or conflicting calibration:" : map prettyPrintCalibration xs
    mapM_ (putStr . render) calsDupl
    error "loadCalibrations: Redundant and/or conflicting calibrations (see above)."
  return calsAll

-- | Calibrate height of a node with given path using the uniform distribution.
--
-- If the upper bound is not given, no upper bound is used.
--
-- For reasons of computational efficiency, the path is not checked for
-- validity. Please do so beforehand using 'calibration'.
--
-- Call 'error' if the path is invalid.
calibrateHard ::
  HasHeight a =>
  Calibration ->
  PriorFunction (Tree e a)
calibrateHard c t
  | h <= a' = 0
  | h >* b = 0
  | otherwise = 1
  where
    a' = fromNonNegative a
    h = fromHeight $ getHeightFromNode p t
    (Interval a b) = calibrationInterval c
    p = calibrationNode c

-- | Calibrate height of a node with given path.
--
-- When the node is in the given bounds, a uniform distribution is used.
--
-- When the node is out of bounds, a one-sided normal distribution with given
-- standard deviation is used. The normal distribution is normalized such that
-- the complete distribution of the constraint is continuous. Use of the normal
-- distribution also ensures that the first derivative is continuous.
--
-- If the upper bound is not given, no upper bound is used.
--
-- For reasons of computational efficiency, the path is not checked for
-- validity. Please do so beforehand using 'calibration'.
--
-- Call 'error' if the path is invalid.
calibrateSoft ::
  HasHeight a =>
  StandardDeviation ->
  Calibration ->
  PriorFunction (Tree e a)
calibrateSoft s c t
  | h <= a' = Exp $ logDensity d (a' - h) - logDensity d 0
  | h >* b = case b of
    Infinity -> 1.0
    Positive b' -> Exp $ logDensity d (h - b') - logDensity d 0
  | otherwise = 1
  where
    a' = fromNonNegative a
    h = fromHeight $ getHeightFromNode p t
    d = normalDistr 0 s
    (Interval a b) = calibrationInterval c
    p = calibrationNode c

-- XXX: Here, we may have to extract the heights first and then check them. Or
-- go through all nodes and check if there is a calibration.

-- | Calibrate nodes of a tree using 'calibrateSoft'.
--
-- Calculate the calibration prior for a given vector of calibrations, the
-- absolute height of the tree, and the tree with relative heights.
--
-- Calibrations can be created using 'calibration' or 'loadCalibrations'. The
-- reason is that finding the nodes on the tree is a slow process not to be
-- repeated at each proposal.
--
-- Call 'error' if:
--
-- - A path is invalid.
--
-- - The height multiplier is zero or negative.
calibrate ::
  HasHeight a =>
  -- | Standard deviation of the calibrations before scaling with the height
  -- multiplier.
  StandardDeviation ->
  V.Vector Calibration ->
  -- | Height multiplier. Useful if working on normalized trees.
  Double ->
  PriorFunction (Tree e a)
calibrate sd cs h t
  | h <= 0 = error "calibrate: Height multiplier is zero or negative."
  | otherwise = V.product $ V.map f cs
  where
    f (Calibration n x i) =
      let i' = if h == 1.0 then i else transformInterval (recip h) i
       in calibrateSoft sd (Calibration n x i') t
