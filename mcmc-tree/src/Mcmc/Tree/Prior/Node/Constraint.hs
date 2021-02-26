{-# LANGUAGE DeriveGeneric #-}

-- |
-- Module      :  Mcmc.Tree.Prior.Node.Constraint
-- Description :  Relative node order constraints and node calibrations
-- Copyright   :  (c) Dominik Schrempf, 2020
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

-- Check if a constraint is valid.
--
-- Return 'Left' if:
--
-- - The younger node is a direct ancestor of the old node.
--
-- - The older node is a direct ancestor of the young node.
validateConstraint :: Constraint -> Either String Constraint
validateConstraint c
  | y `isPrefixOf` o = Left $ getErrMsg "Younger node is direct ancestor of older node (?)."
  | o `isPrefixOf` y = Left $ getErrMsg "No need to constrain older node which is direct ancestor of younger node."
  | otherwise = Right c
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

-- | Load constraints from file.
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
loadConstraints :: Tree e Name -> FilePath -> IO (V.Vector Constraint)
loadConstraints t f = do
  d <- BL.readFile f
  let mr = decode NoHeader d :: Either String (V.Vector ConstraintData)
      cds = either error id mr
  return $ V.map (constraintDataToConstraint t) cds

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
