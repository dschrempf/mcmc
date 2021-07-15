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
    calibrateSoftF,
    calibrate,
  )
where

import Control.Lens
import Control.Monad
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Csv hiding (Name)
import Data.Function
import Data.List
import qualified Data.Vector as V
import ELynx.Tree hiding (partition)
import GHC.Generics
import Mcmc.Prior hiding (positive)
import Mcmc.Statistics.Types
import Mcmc.Tree.Lens
import Mcmc.Tree.Mrca
import Mcmc.Tree.Types
import Text.Read

-- | Non-negative number.
newtype NonNegative a = NonNegative {fromNonNegative :: a}
  deriving (Eq)

nonNegative :: (Ord a, Num a) => a -> NonNegative a
nonNegative x
  | x < 0 = error "nonNegative: Negative value."
  | otherwise = NonNegative x

instance (RealFloat a, Read a) => Read (NonNegative a) where
  readPrec = nonNegative <$> readPrec

instance Show a => Show (NonNegative a) where
  showsPrec p (NonNegative x) = showsPrec p x

-- | Positive number or infinity.
data ExtendedPositive a = Positive a | Infinity
  deriving (Eq)

positive :: (Ord a, Num a) => a -> ExtendedPositive a
positive x
  | x <= 0 = error "positive: Zero or negative value."
  | otherwise = Positive x

positiveReadPrec :: (Ord a, Num a, Read a) => ReadPrec (ExtendedPositive a)
positiveReadPrec = positive <$> readPrec

infinityReadPrec :: ReadPrec (ExtendedPositive a)
infinityReadPrec = do
  Ident "Infinity" <- lexP
  return Infinity

instance (Ord a, Num a, Read a) => Read (ExtendedPositive a) where
  readPrec = positiveReadPrec <++ infinityReadPrec

instance Show a => Show (ExtendedPositive a) where
  showsPrec p (Positive x) = showsPrec p x
  showsPrec _ Infinity = showString "Infinity"

-- | Open interval \((a,b)\) with \(a < b\), \(a \in [0, \infty)\) and \(b \in
-- (0, \infty]\).
data Interval a = Interval (NonNegative a) (ExtendedPositive a)
  deriving (Eq)

instance Show a => Show (Interval a) where
  show (Interval a b) = "(" ++ show a ++ ", " ++ show b ++ ")"

-- | Specify a lower and an upper bound.
properInterval :: (Ord a, Num a) => LowerBoundary a -> UpperBoundary a -> Interval a
properInterval a b
  | a < b = Interval (nonNegative a) (positive b)
  | otherwise = error "properInterval: Left bound equal or larger right bound."
{-# SPECIALIZE properInterval :: Double -> Double -> Interval Double #-}

-- | Specify a lower bound only. The upper bound is set to infinity.
lowerBoundOnly :: (Ord a, Num a) => LowerBoundary a -> Interval a
lowerBoundOnly a = Interval (nonNegative a) Infinity
{-# SPECIALIZE lowerBoundOnly :: Double -> Interval Double #-}

-- | Transform an interval by applying a multiplicative change.
--
-- Useful when the tree is normalized and height values have to be converted
-- from relative heights to absolute heights.
transformInterval :: RealFloat a => a -> Interval a -> Interval a
transformInterval x (Interval a b)
  | x <= 0 = error "transform: Multiplier is zero or negative."
  | otherwise = Interval a' b'
  where
    a' = NonNegative $ x * fromNonNegative a
    b' = case b of
      Positive bVal -> Positive $ x * bVal
      Infinity -> Infinity
{-# SPECIALIZE transformInterval :: Double -> Interval Double -> Interval Double #-}

-- No number is bigger than a non-existing upper bound..
(>*) :: (Ord a, Fractional a) => a -> ExtendedPositive a -> Bool
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
data Calibration a = Calibration
  { calibrationName :: String,
    calibrationNodePath :: Path,
    calibrationNodeIndex :: Int,
    calibrationInterval :: Interval a
  }
  deriving (Eq, Show)

prettyPrintCalibration :: Show a => Calibration a -> String
prettyPrintCalibration (Calibration n p i l) =
  "Calibration: "
    <> n
    <> " with path "
    <> show p
    <> ", index "
    <> show i
    <> ", and interval "
    <> show l
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
  Interval b ->
  Calibration b
calibration t n xs = Calibration n p i
  where
    err msg = error $ "calibration: " ++ n ++ ": " ++ msg
    p = either err id $ mrca xs t
    -- NOTE: Identifying the tree multiple times may be slow when creating many
    -- calibrations. But this is only done once in the beginning.
    i = label $ getSubTreeUnsafe p $ identify t
{-# SPECIALIZE calibration ::
  (Ord a, Show a) =>
  Tree e a ->
  String ->
  [a] ->
  Interval Double ->
  Calibration Double
  #-}

-- Used to decode the CSV file.
data CalibrationData a
  = CalibrationData
      String -- Calibration name.
      String -- Leaf a.
      String -- Leaf b.
      a -- Lef boundary.
      (Maybe a) -- Maybe right boundary.
  deriving (Generic, Show)

instance FromField a => FromRecord (CalibrationData a)

calibrationDataToCalibration :: (Ord a, Num a) => Tree e Name -> CalibrationData a -> Calibration a
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
loadCalibrations :: Tree e Name -> FilePath -> IO (V.Vector (Calibration Double))
loadCalibrations t f = do
  d <- BL.readFile f
  let mr = decode NoHeader d :: Either String (V.Vector (CalibrationData Double))
      cds = either error id mr
  when (V.null cds) $ error $ "loadCalibrations: No calibrations found in file: " <> f <> "."
  let calsAll = V.map (calibrationDataToCalibration t) cds
  -- Check for duplicates and conflicts.
  let calsDupl = findDupsBy ((==) `on` calibrationNodePath) $ V.toList calsAll
  if null calsDupl
    then putStrLn "No duplicates and no conflicting calibrations have been detected."
    else do
      -- Calibrations could also be removed. But then, which one should be removed?
      let render xs =
            unlines $
              "Redundant and/or conflicting calibration:" : map prettyPrintCalibration xs
      mapM_ (putStr . render) calsDupl
      error "loadCalibrations: Duplicates and/or conflicting calibrations have been detected."
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
  RealFloat a =>
  Calibration a ->
  PriorFunctionG (HeightTree a) a
calibrateHard c (HeightTree t)
  | h <= a' = 0.0
  | h >* b = 0.0
  | otherwise = 1.0
  where
    a' = realToFrac $ fromNonNegative a
    h = t ^. subTreeAtL p . branchL
    (Interval a b) = calibrationInterval c
    p = calibrationNodePath c

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
  RealFloat a =>
  StandardDeviation a ->
  Calibration a ->
  PriorFunctionG (HeightTree a) a
calibrateSoft s c (HeightTree t) = calibrateSoftF s l h
  where
    p = calibrationNodePath c
    h = t ^. subTreeAtL p . branchL
    l = calibrationInterval c

-- | See 'calibrateSoft'.
calibrateSoftF :: RealFloat a => StandardDeviation a -> Interval a -> PriorFunctionG a a
calibrateSoftF s (Interval a' b) h
  | h <= a = d (a - h) / d 0
  | h >* b = case b of
    Infinity -> 1.0
    Positive b' -> d (h - b') / d 0
  | otherwise = 1
  where
    a = fromNonNegative a'
    d = normal 0 s

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
  RealFloat a =>
  -- | Standard deviation of the calibrations before scaling with the height
  -- multiplier.
  --
  -- NOTE: The same standard deviation is used for all calibrations.
  StandardDeviation a ->
  V.Vector (Calibration a) ->
  -- | Height multiplier of tree. Useful when working on normalized trees.
  a ->
  PriorFunctionG (HeightTree a) a
calibrate sd cs h t
  | h <= 0 = error "calibrate: Height multiplier is zero or negative."
  | otherwise = V.product $ V.map f cs
  where
    f (Calibration n x i l) =
      let l' = if h == 1.0 then l else transformInterval (recip h) l
       in calibrateSoft sd (Calibration n x i l') t
