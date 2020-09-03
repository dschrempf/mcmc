-- |
-- Module      :  Mcmc.Tree.Prior.MolecularClock
-- Description :  Molecular clock priors
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Sep  3 14:20:55 2020.
module Mcmc.Tree.Prior.MolecularClock
  ( uncorrelatedGamma,
  )
where

import ELynx.Data.Tree
import Mcmc.Prior
import Numeric.Log

-- | Uncorrelated gamma prior.
--
-- The first tree stores the times. The second tree stores the rates.
--
-- For a given time range @dt@, the rate is distributed according to a gamma
-- distribution with mean @dt@ and variance @dt^2@.
--
-- The root branch is ignored!
uncorrelatedGamma :: Tree Double a -> Tree Double b -> Log Double
uncorrelatedGamma t r = product' [gamma m v r | (m, v, r) <- zip3 ts vs rs]
  where
    ts = tail $ branches t
    vs = map (** 2) ts
    rs = tail $ branches r
