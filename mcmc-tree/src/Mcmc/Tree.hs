-- |
-- Module      :  Mcmc.Tree
-- Description :  Markov chain Monte Carlo sampling on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri Jul  3 10:09:17 2020.
module Mcmc.Tree
  ( Tree,
    getLens,
    newick,
    oneNewick,
    manyNewick,
  )
where

import Mcmc.Tree.Tree
import Mcmc.Tree.Newick
