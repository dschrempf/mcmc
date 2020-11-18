-- |
-- Module      :  Mcmc.Tree.Proposal.Common
-- Description :  Common functions used by all tree proposals
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Nov  4 11:53:16 2020.
module Mcmc.Tree.Proposal.Common
  ( lengthL,
    truncatedNormalSample,
  )
where

import Control.Lens
import ELynx.Tree
import Mcmc.Tree.Lens
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.TruncatedNormal
import System.Random.MWC

-- | Which lens to use?
--
-- Used for debugging, `lengthE` is around ten percent slower than `lengthU`.
lengthL :: Lens' Length Double
-- lengthL = lengthE
lengthL = lengthU

-- | A very specific function that samples a delta value from the truncated
-- normal distribution with given bounds [a,b] and also computes the required
-- factor of the Metropolis-Hastings proposal ratio.
--
-- NO JACOBIAN IS COMPUTED, because we do not know how the proposal will be used.
truncatedNormalSample ::
  -- Standard deviation.
  Double ->
  -- Tuning parameter.
  Double ->
  -- Left bound.
  Double ->
  -- Right bound.
  Double ->
  GenIO ->
  IO (Double, Log Double)
truncatedNormalSample s t a b g = do
  let s' = t * s
      d = either error id $ truncatedNormalDistr 0 s' a b
  u <- genContinuous d g
  -- Compute Metropolis-Hastings factor.
  let a' = a - u
      b' = b - u
      d' = either error id $ truncatedNormalDistr 0 s' a' b'
      qXY = Exp $ logDensity d u
      qYX = Exp $ logDensity d' (- u)
  -- NO JACOBIAN IS COMPUTED.
  return (u, qYX / qXY)
