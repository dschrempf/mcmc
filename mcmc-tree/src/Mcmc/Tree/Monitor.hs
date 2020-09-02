-- |
-- Module      :  Mcmc.Tree.Monitor
-- Description :  Monitors for trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Sat Jul 25 14:06:14 2020.
module Mcmc.Tree.Monitor
  ( monitorTree,
  )
where

import ELynx.Data.Tree
import ELynx.Export.Tree.Newick
import Mcmc

-- | Monitor a tree in Newick format.
monitorTree ::
  (Measurable e, Named a) =>
  -- | Name.
  String ->
  MonitorParameter (Tree e a)
monitorTree n =
  MonitorParameter
    n
    (toNewickBuilder . measurableToPhyloTree)
