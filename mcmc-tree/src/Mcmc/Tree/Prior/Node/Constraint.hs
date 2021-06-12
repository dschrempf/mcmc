{-# LANGUAGE DeriveGeneric #-}

-- |
-- Module      :  Mcmc.Tree.Prior.Node.Constraint
-- Description :  Relative node order constraints and node calibrations
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jul 27 10:49:11 2020.
module Mcmc.Tree.Prior.Node.Constraint
  ( -- * Constraints
    Constraint (..),
    constraint,
    loadConstraints,
    constrainHardUnsafe,
    constrainSoftUnsafe,
    constrainUnsafe,
  )
where

import Control.Monad
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Csv hiding (Name)
import Data.List
import qualified Data.Vector as V
import ELynx.Tree
import GHC.Generics
import Mcmc.Chain.Chain
import Mcmc.Statistics.Types
import Mcmc.Tree.Mrca
import Mcmc.Tree.Prior.Node.Common
import Mcmc.Tree.Types
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal

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

-- Is the left node an ancestor of the right node?
isAncestor :: Eq a => [a] -> [a] -> Bool
isAncestor = isPrefixOf

-- Is the left node a descdant of the right node?
isDescendant :: Eq a => [a] -> [a] -> Bool
isDescendant = flip isPrefixOf

-- Relationship of two nodes.
data Relationship
  = Equal
  | -- Same as RightIsDescendentOfLeft.
    LeftIsAncestorOfRight
  | -- Same as RightIsAncestorOfLeft.
    LeftIsDescendantOfRight
  | Unrelated
  deriving (Eq, Show, Read)

-- Fast function avoiding two consecutive uses of `isPrefixOf`.
--
-- Are the two nodes direct descent of each other? Basically check both:
-- 'isAncestor' and 'isDescendent'.
areDirectDescendants :: Eq a => [a] -> [a] -> Relationship
areDirectDescendants [] [] = Equal
areDirectDescendants (x : xs) (y : ys)
  | x == y = areDirectDescendants xs ys
  | otherwise = Unrelated
areDirectDescendants (_ : _) [] = LeftIsDescendantOfRight
areDirectDescendants [] (_ : _) = LeftIsAncestorOfRight

-- Check if a constraint is valid.
--
-- See also 'validateConstraints' and 'validateConstraintVector' which perform
-- checks on a multiple possibly redundant or conflicting constraints.
--
-- Return 'Left' if:
--
-- - The younger node is a direct ancestor of the old node.
--
-- - The older node is a direct ancestor of the young node.
--
-- NOTE: Any node is a direct ancestor of itself, and so, bogus constraints
-- including the same node twice are also filtered out.
validateConstraint :: Constraint -> Either String Constraint
validateConstraint c = case areDirectDescendants y o of
  Equal ->
    Left $
      getErrMsg "Bogus constraint; both nodes are equal (?)."
  LeftIsAncestorOfRight ->
    Left $
      getErrMsg "Bogus constraint; younger node is direct ancestor of older node (?)."
  LeftIsDescendantOfRight ->
    Left $
      getErrMsg "Redundant constraint; old node is direct ancestors of young node."
  Unrelated -> Right c
  where
    n = constraintName c
    y = constraintYoungNode c
    o = constraintOldNode c
    getErrMsg msg = "validateConstraint: " ++ show n ++ ": " ++ msg

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
    err msg = error $ "constraint: " ++ show n ++ ": " ++ msg
    y = either err id $ mrca ys t
    o = either err id $ mrca os t

-- XXX: Maybe use a vector (but we may just fold over it then a list is fine).
--
-- type Constraints = [Constraint]

data ConstraintData = ConstraintData String String String String String
  deriving (Generic, Show)

instance FromRecord ConstraintData

constraintDataToConstraint :: Tree e Name -> ConstraintData -> Constraint
constraintDataToConstraint t (ConstraintData n yL yR oL oR) =
  constraint t n [f yL, f yR] [f oL, f oR]
  where
    f = Name . BL.pack

-- Given a left constraint, another right constraint can be:
data Property
  = RightIsRedundant Constraint Constraint
  | RightIsConflicting Constraint Constraint
  | RightIsFine Constraint Constraint
  deriving (Eq, Show, Read)

isFine :: Property -> Bool
isFine (RightIsFine _ _) = True
isFine _ = False

describeProp :: Property -> String
describeProp (RightIsRedundant l r) =
  "Constraint " <> constraintName r <> " is redundant given constraint " <> constraintName l <> "."
describeProp (RightIsConflicting l r) =
  "Constraint " <> constraintName r <> " conflicts with constraint " <> constraintName l <> "."
describeProp (RightIsFine l r) =
  "Constraints " <> constraintName l <> " and " <> constraintName r <> " are fine."

-- Given the left constraint, check if the right constraint is redundant or
-- conflicting.
--
-- See also 'validateConstraint' which performs checks on single constraints.
--
-- This validation includes a complicated list of tests checking for
-- redundancies or conflicts between two constraints.
--
-- Let D(x,y) be true iff x is a descendant of y (see 'isDescendent').
--
-- Let A(x,y) be true iff x is an ancestor of y (see 'isAncestor').
--
-- Given the left constraint a < b, the right constraint c < d is
--
-- - redundant iff: D(c,a) AND A(d,b)
--
-- - conflicting iff: A(c,b) AND D(d,a)
--
-- - otherwise fine.
--
-- Did I miss any other tests?
--
validateConstraints :: Constraint -> Constraint -> Property
validateConstraints l@(Constraint _ a b) r@(Constraint _ c d)
  | (c `isDescendant` a) && (d `isAncestor` b) = RightIsRedundant l r
  | (c `isAncestor` b) && (d `isDescendant` a) = RightIsConflicting l r
  | otherwise = RightIsFine l r

-- Pairwise comparison with 'validateConstraints'.
validateConstraintVector :: V.Vector Constraint -> V.Vector Property
validateConstraintVector xs = do
  l <- xs
  r <- xs
  guard (l /= r)
  return $ validateConstraints l r

-- | Load and validate constraints from file.
--
-- The constraint file is a comma separated values (CSV) file with rows of the
-- following format:
--
-- @
-- ConstraintName,YoungLeafA,YoungLeafB,OldLeafA,OldLeafB
-- @
--
-- The young and old nodes are uniquely defined as the most recent common
-- ancestors (MRCA) @YoungLeafA@ and @YoungLeafB@, as well as @OldLeafA@ and
-- @OldLeafB@.
--
-- The following line defines a constraint where the ancestor of leaves A and B
-- is younger than the ancestor of leaves C and D:
--
-- @
-- ExampleConstraint,A,B,C,D
-- @
--
-- Call 'error' if
--
-- - The file contains errors.
--
-- - An MRCA cannot be found.
--
-- - Redundant or conflicting constraints are found.
loadConstraints :: Tree e Name -> FilePath -> IO (V.Vector Constraint)
loadConstraints t f = do
  d <- BL.readFile f
  let mr = decode NoHeader d :: Either String (V.Vector ConstraintData)
      cds = either error id mr
  when (V.null cds) $ error $ "loadConstraints: No constraints found in file: " <> f <> "."
  let cons = V.map (constraintDataToConstraint t) cds
  -- Check for conflicting constraints.
  let badCons = V.filter (not . isFine) $ validateConstraintVector cons
  if V.null badCons
    then return cons
    else do
      V.sequence_ $ V.map (putStrLn . describeProp) badCons
      error "loadConstraints: Constraints are erroneous (see above)."

-- | Hard constrain order of nodes with given paths.
--
-- A truncated, improper uniform distribution is used.
--
-- For reasons of computational efficiency, the paths are not checked for
-- validity. Please do so beforehand using 'constraint'.
constrainHardUnsafe ::
  HasHeight a =>
  Constraint ->
  PriorFunction (Tree e a)
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
-- validity. Please do so beforehand using 'constraint'.
constrainSoftUnsafe ::
  HasHeight a =>
  StandardDeviation ->
  Constraint ->
  PriorFunction (Tree e a)
constrainSoftUnsafe s c t
  | hY < hO = 1
  | otherwise = Exp $ logDensity d (hY - hO) - logDensity d 0
  where
    hY = fromHeight $ getHeightFromNodeUnsafe y t
    hO = fromHeight $ getHeightFromNodeUnsafe o t
    d = normalDistr 0 s
    y = constraintYoungNode c
    o = constraintOldNode c

-- XXX: Here, we may have to extract the heights first and then check them. Or
-- go through all nodes and check if there is a calibration.

-- | Constrain nodes of a tree using 'constrainSoftUnsafe'.
--
-- Calculate the constraint prior for a given vector of constraints, and a
-- tree with relative heights.
--
-- Constraints can be created using 'constraint' or 'loadConstraints'. The
-- reason is that finding the nodes on the tree is a slow process not to be
-- repeated at each proposal.
--
-- Call 'error' if a path is invalid.
constrainUnsafe ::
  HasHeight a =>
  StandardDeviation ->
  V.Vector Constraint ->
  PriorFunction (Tree e a)
constrainUnsafe sd cs t = V.product $ V.map f cs
  where
    f x = constrainSoftUnsafe sd x t
