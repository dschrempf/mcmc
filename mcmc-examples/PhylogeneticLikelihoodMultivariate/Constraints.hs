{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Constraints
-- Description :  Constrain node order
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Oct 22 17:55:06 2020.
module Constraints
  ( Constraint,
    constraints,
    getConstraints,
  )
where

import ELynx.Tree
import Mcmc.Tree
import Numeric.Log

-- | Constraints define node orders.
--
-- For example, @("MyConstraint", YOUNG, OLD)@ ensures node @YOUNG@ to be
-- younger than node @OLD@, and gives the constraint the name @MyConstraint@.
type Constraint = (String, Path, Path)

-- | Constraint prior.
--
-- For a given set of constraints and the relative time tree, calculate the
-- constraint prior.
--
-- The constraints have to be pre-computed with 'getConstraints'. The reason is
-- that finding the nodes on the tree is a slow process that should not be
-- repeated.
constraints :: HasHeight a => [Constraint] -> Tree e a -> [Log Double]
constraints xs t =
  [constrainSoft 1e-3 y o t | (_, y, o) <- xs]

-- | Find and constrain the constrained nodes on the tree.
getConstraints :: Tree e Name -> [Constraint]
getConstraints t = [("hornwort->fern", young, old)]
  where
    young = mrcaUnsafe ["Lindsaea_linearis", "Polystichum_acrostichoides"] t
    old = mrcaUnsafe ["Phaeomegaceros_coriaceus", "Nothoceros_aenigmaticus"] t
