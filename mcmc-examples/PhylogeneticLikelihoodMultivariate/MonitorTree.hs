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
  )
where

import Data.Bifunctor
import qualified Data.ByteString.Conversion.From as L
import qualified Data.ByteString.Lazy.Char8 as L
import Data.Maybe
import qualified Data.Text.Lazy.Builder as B
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
  MonitorParameter n (B.fromText . toText . toNewick . lengthToPhyloTree . first Length)
  where
    toText s = fromMaybe (error "conversion failed") $ L.fromByteString $ L.toStrict s
