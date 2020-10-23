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
import Data.Maybe
import ELynx.Tree
import Mcmc.Tree
import Numeric.Log

-- | Constraints define node orders.
--
-- (a, b) ensures node a to be younger than node b.
type Constraint = (Path, Path)

-- | Constraint prior.
--
-- TODO: Write documentation;
constraints :: [Constraint] -> Tree Double Double -> [Log Double]
constraints xs t =
  [constrainSoft 1e-4 y o t | (y, o) <- xs]

-- | The constraint induced by a horizontal gene transfer does not contradict
-- the node order of the substitution-like tree obtained from the sequences.
getConstraints :: Tree e BS.ByteString -> [Constraint]
getConstraints t = [(young, old)]
  where
    young =
      fromMaybe
        (error "constrainedNodes: No MRCA young.")
        -- (mrca [213, 200])
        (mrca ["Pt_vittata", "Po_acrosti"] t)
    old =
      fromMaybe
        (error "constrainedNodes: No MRCA old.")
        -- (mrca [144, 143, 142])
        (mrca ["Me_tosanus", "Me_vincent", "No_aenigma"] t)
