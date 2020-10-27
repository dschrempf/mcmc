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

import qualified Data.ByteString.Char8 as BS
import ELynx.Tree
import Mcmc.Tree
import Numeric.Log

-- | Constraints define node orders.
--
-- @("name", a, b)@ ensures node @a@ to be younger than node @b@, and gives the
-- constraint the name @name@.
type Constraint = (String, Path, Path)

-- | Constraint prior.
--
-- For a given set of constraints and the relative time tree, calculate the
-- constraint prior.
--
-- The constraints have to be pre-computed with 'getConstraints'. The reason is
-- that finding the nodes on the tree is a slow process that should not be
-- repeated.
constraints :: [Constraint] -> Tree Double Double -> [Log Double]
constraints xs t =
  [constrainSoft 1e-3 y o t | (_, y, o) <- xs]

-- | Find and constrain the constrained nodes on the tree.
getConstraints :: Tree e BS.ByteString -> [Constraint]
getConstraints t = [("hornwort->fern", young, old)]
  where
    young =
      fromMaybe
        (error "constrainedNodes: No MRCA young.")
        (mrca ["Lindsaea_linearis", "Polystichum_acrostichoides"] t)
    old =
      fromMaybe
        (error "constrainedNodes: No MRCA old.")
        (mrca ["Phaeomegaceros_coriaceus", "Nothoceros_aenigmaticus"] t)
