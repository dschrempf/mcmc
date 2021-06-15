-- |
-- Module      :  Mcmc.Tree.Prior.Node.Common
-- Description :  Commonly used functions of tree node priors
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Thu Feb 25 15:34:52 2021.
module Mcmc.Tree.Prior.Node.Common
  ( getHeightFromNode,
  )
where

import Control.Lens
import ELynx.Tree
import Mcmc.Tree.Lens
import Mcmc.Tree.Types

-- | Get the height of the node at path on the tree.
--
-- Call 'error' if the path is invalid.
getHeightFromNode :: HasHeight a => Path -> Tree e a -> Height
getHeightFromNode p t = t ^. subTreeAtL p . labelL . hasHeightL
