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
  ( truncatedNormalSample,
  )
where

import Control.Monad
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.TruncatedNormal
import System.Random.MWC

-- | A very specific function that samples a delta value from the truncated
-- normal distribution with given bounds [a,b] and also computes the required
-- factor of the Metropolis-Hastings-Green (MHG) proposal ratio.
--
-- NO JACOBIAN IS COMPUTED, because we do not know how the proposal will be used.
truncatedNormalSample ::
  -- | Mean.
  Double ->
  -- | Standard deviation.
  Double ->
  -- | Tuning parameter.
  Double ->
  -- | Left bound.
  Double ->
  -- | Right bound.
  Double ->
  GenIO ->
  -- | (NewValue, MHGRatioWithoutJacobian)
  IO (Double, Log Double)
truncatedNormalSample m s t a b g = do
  let s' = t * s
      d = truncatedNormalDistr m s' a b
  u <- genContinuous d g
  let err = error $
        "truncatedNormalSample: Value "
          <> show u
          <> " out of bounds ["
          <> show a
          <> ","
          <> show b
          <> "]."
  when (a > u || b < u) err
  -- Compute Metropolis-Hastings-Green factor.
  let d' = truncatedNormalDistr u s' a b
      qXY = Exp $ logDensity d u
      qYX = Exp $ logDensity d' m
  -- NO JACOBIAN IS COMPUTED.
  return (u, qYX / qXY)
