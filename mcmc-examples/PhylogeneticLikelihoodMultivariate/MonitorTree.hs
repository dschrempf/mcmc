-- |
-- Module      :  MonitorTree
-- Description :  Monitors for trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Sat Jul 25 14:06:14 2020.
module MonitorTree
  ( monitorTree,
    monitorTreeWith,
  )
where

import Data.Bifunctor
import ELynx.Data.Tree
import ELynx.Export.Tree.Newick
import Mcmc

-- | Monitor a tree in Newick format.
monitorTree ::
  Named a =>
  -- | Name.
  String ->
  MonitorParameter (Tree Double a)
monitorTree n =
  MonitorParameter
    n
    (toNewickBuilder . lengthToPhyloTree . first Length)

-- | Monitor a tree in Newick format.
--
-- Apply the given function before executing the monitor.
monitorTreeWith ::
  (Measurable f, Named b) =>
  -- | Function to apply to tree.
  (Tree e a -> Tree f b) ->
  -- | Name.
  String ->
  MonitorParameter (Tree e a)
monitorTreeWith f n =
  MonitorParameter
    n
    (toNewickBuilder . lengthToPhyloTree . first (Length . getLen) . f)
    -- TODO: This construct is weird. Probably provide something like
    --
    -- fromMeasurable :: Measurable e -> e -> Phylo
