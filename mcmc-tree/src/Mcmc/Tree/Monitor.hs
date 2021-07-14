-- |
-- Module      :  Mcmc.Tree.Monitor
-- Description :  Monitors for trees
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Sat Jul 25 14:06:14 2020.
module Mcmc.Tree.Monitor
  ( monitorTree,
    monitorTreeWith,
  )
where

import Data.Bifunctor
import Data.Maybe
import ELynx.Tree
import Mcmc

-- | Monitor a tree in Newick format.
monitorTree ::
  (HasLength e, HasName a) =>
  -- | Name.
  String ->
  MonitorParameter (Tree e a)
monitorTree n = MonitorParameter n (toNewickBuilder . lengthToPhyloTree)

-- | Monitor a tree in Newick format.
monitorTreeWith ::
  -- | Node names in pre-order.
  [Name] ->
  -- | Name.
  String ->
  MonitorParameter (Tree Double ())
monitorTreeWith ns n =
  MonitorParameter
    n
    (toNewickBuilder . lengthToPhyloTree . checkLengths . setNames)
  where
    checkLengths = first (either error id . toLength)
    setNames =
      fromMaybe (error "monitorTreeWith: Too few or too many node labels.") . setLabels ns
