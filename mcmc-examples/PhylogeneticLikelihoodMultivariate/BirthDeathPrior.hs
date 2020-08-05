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
  ( birthDeathWith,
  )
where

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
birthDeathWith ::
  -- | Birth rate.
  Double ->
  -- | Death rate.
  Double ->
  Tree Double a ->
  Log Double
-- TODO: Check for critical and nearly critical process.
birthDeathWith l m
  | l < 0.0 = error "birthDeathWith: Birth rate lambda is negative."
  | m < 0.0 = error "birthDeathWith: Death rate mu is negative."
  | epsNearCritical > abs (l - m) = error "birthDeathWith: Birth and death rate are too close."
  | otherwise = Exp . log . fst . birthDeathWith' l m

birthDeathWith' :: Double -> Double -> Tree Double a -> (Double, Double)
birthDeathWith' la mu (Node br _ []) = computeDE la mu br 0
birthDeathWith' la mu (Node br _ [l, r]) = (dT * dL * dR * la, eT)
  where
    (dL, eL) = birthDeathWith' la mu l
    (dR, eR) = birthDeathWith' la mu r
    (dT, eT) = computeDE la mu br (eL * eR)
birthDeathWith' _ _ _ = error "birthDeathWith: Tree is not bifurcating."

-- Tests

-- testTree1 :: Tree Double ()
-- testTree1 = Node 1.0 () []

-- >>> birthDeathWith 1.2 3.2 testTree1
-- 5.8669248906043234e-2

-- testTree2 :: Tree Double ()
-- testTree2 = Node 0.0 () [Node 0.4 () [], Node 0.2 () [Node 0.2 () [], Node 0.2 () []]]

-- >>> birthDeathWith 1.2 3.2 testTree2
-- 3.978845396350806e-2

-- XXX: There are differences in the conditions. The point process conditions on
-- the time of origin, and on the number of leaves. The dynamic programming
-- approach used above only conditions on the time of origin.

-- -- TODO.
-- -- Assume that node labels denote node heights.
-- birthDeathPointProcessWith ::
--   -- | Birth rate.
--   Double ->
--   -- | Death rate.
--   Double ->
--   Tree Double Double ->
--   Log Double
-- -- TODO: Check for critical and nearly critical process.
-- birthDeathPointProcessWith l m t
--   | l < 0.0 = error "birthDeathWith: Birth rate lambda is negative."
--   | m < 0.0 = error "birthDeathWith: Death rate mu is negative."
--   | otherwise = Exp $ logDensity d 0
--   where
--     d = BDD (label t + branch t) l m
