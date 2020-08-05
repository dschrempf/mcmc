-- |
-- Module      :  BirthDeathPrior
-- Description :  Calculate probability of a tree assuming birth and death process
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Aug  4 20:37:11 2020.
module BirthDeathPrior
  ( birthDeath,
  )
where

import Data.Bifunctor
import ELynx.Data.Tree
import Numeric.Log

-- Compute probabilities D and E at the top of the branch.
--
-- Correct results:
-- >>> computeDE 1.2 3.2 1.0 0.3
-- (7.283127121752474e-2,0.9305035687810801)
computeDE ::
  -- Lambda.
  Double ->
  -- Mu.
  Double ->
  -- Branch length.
  Double ->
  -- E at bottom of branch.
  Double ->
  -- D, E.
  (Double, Double)
computeDE la mu dt e0 = (a / b / b, c / b)
  where
    d = la - mu
    x = exp (- d * dt)
    y = (mu - e0 * la) * x
    e' = e0 - 1.0
    laE' = la * e'
    a = d * d * x
    b = laE' + y
    c = mu * e' + y
{-# INLINE computeDE #-}

-- Require near critical process if birth and death rates are closer than this value.
epsNearCritical :: Double
epsNearCritical = 1e-5

-- -- Require critical process if birth and death rates are closer than this value.
-- eps :: Double
-- eps = 1e-12

-- | Birth and death prior.
--
-- The sampling rate is 1.0; i.e., the extinction probability of leaves is 0.0.
--
-- XXX: This prior does not condition on survival.
--
-- XXX: This prior conditions on the origin, and not on the root node. This
-- affects: (1) the time, and (2) the split at the root, because the prior has
-- an additional multiplicative factor.
--
-- The aforementioned points are just multiplicative factors that don't
-- influence the stationary distribution of the MCMC run. However, they do
-- change the prior, and hence, the posterior.
birthDeath ::
  -- | Birth rate.
  Double ->
  -- | Death rate.
  Double ->
  Tree Double a ->
  Log Double
-- TODO: Check for critical and nearly critical process.
birthDeath la mu
  | la < 0.0 = error "birthDeath: Birth rate lambda is negative."
  | mu < 0.0 = error "birthDeath: Death rate mu is negative."
  | epsNearCritical > abs (la - mu) = error "birthDeath: Birth and death rate are too close."
  | otherwise = fst . birthDeath' la mu (Exp $ log la)

birthDeath' :: Double -> Double -> Log Double -> Tree Double a -> (Log Double, Double)
birthDeath' la mu _ (Node br _ []) = first (Exp . log) $ computeDE la mu br 0
birthDeath' la mu logLa (Node br _ [l, r]) = ((Exp $ log dT) * dL * dR * logLa, eT)
  where
    (dL, eL) = birthDeath' la mu logLa l
    (dR, eR) = birthDeath' la mu logLa r
    (dT, eT) = computeDE la mu br (eL * eR)
birthDeath' _ _ _ _ = error "birthDeath: Tree is not bifurcating."

-- Tests

-- testTree1 :: Tree Double ()
-- testTree1 = Node 1.0 () []

-- >>> birthDeath 1.2 3.2 testTree1
-- 5.8669248906043234e-2

-- testTree2 :: Tree Double ()
-- testTree2 = Node 0.0 () [Node 0.4 () [], Node 0.2 () [Node 0.2 () [], Node 0.2 () []]]

-- >>> birthDeath 1.2 3.2 testTree2
-- 3.978845396350806e-2

-- XXX: There are differences in the conditions. The point process conditions on
-- the time of origin, and on the number of leaves. The dynamic programming
-- approach used above only conditions on the time of origin.

-- Question: Is the birth death distribution of the point process adequate, if
-- the topology is given?
--
-- How do we handle the critical, or nearly critical process? Is it fast enough
-- to combine the three distributions into one, and check during calculation of
-- the density? Or is it better to check for the three different possibilities
-- here, and then use the adequate distribution (status quo)?

-- -- TODO.
-- -- Assume that node labels denote node heights.
-- birthDeathPointProcess ::
--   -- | Birth rate.
--   Double ->
--   -- | Death rate.
--   Double ->
--   Tree Double Double ->
--   Log Double
-- -- TODO: Check for critical and nearly critical process.
-- birthDeathPointProcess l m t
--   | l < 0.0 = error "birthDeath: Birth rate lambda is negative."
--   | m < 0.0 = error "birthDeath: Death rate mu is negative."
--   | otherwise = Exp $ logDensity d 0
--   where
--     d = BDD (label t + branch t) l m