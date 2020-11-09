-- |
-- Module      :  Mcmc.Tree.Prior.BirthDeath
-- Description :  Calculate probability of a tree assuming birth and death process
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Aug  4 20:37:11 2020.
module Mcmc.Tree.Prior.BirthDeath
  ( birthDeath,
  )
where

import Data.Bifunctor
-- import Debug.Trace
import ELynx.Tree
import Numeric.Log

-- Compute probabilities D and E at the top of the branch.
--
-- E and D are Eqs [1] and [2] in Stadler, T. Mammalian phylogeny reveals recent
-- diversification rate shifts. Proceedings of the National Academy of Sciences
-- 108, 6187–6192 (2011), respectively.
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
  -- E at bottom of branch. Storage in log domain not necessary.
  Double ->
  -- D, E. Storage in log domain not necessary.
  (Double, Double)
computeDE la mu dt e0 = (a / b / b, c / b)
  where
    d = la - mu
    x = exp (- d * dt)
    y = (mu - e0 * la) * x
    -- Negative survival probability rho.
    --
    -- E0 = 1 - rho
    -- E0 - 1 = -rho
    e' = e0 - 1.0
    laE' = la * e'
    a = d * d * x
    b = laE' + y
    c = mu * e' + y
{-# INLINE computeDE #-}

-- Compute probabilities D and E at the top of the branch for la ~= mu.
--
-- Correct results:
-- >>> computeDENearCritical 1.2 3.2 1.0 0.3
-- (7.283127121752474e-2,0.9305035687810801)
computeDENearCritical ::
  -- Lambda.
  Double ->
  -- Mu.
  Double ->
  -- Branch length.
  Double ->
  -- E at bottom of branch. Storage in log domain not necessary.
  Double ->
  -- D, E. Storage in log domain not necessary.
  (Double, Double)
computeDENearCritical la mu dt e0 = (a / b / b, c / b)
  where
    d = la - mu
    y = mu - e0 * la
    a = 1 - d * dt
    c = e0 + y * dt
    b = 1 + y * dt
{-# INLINE computeDENearCritical #-}

-- Require near critical process if birth and death rates are closer than this value.
epsNearCritical :: Double
epsNearCritical = 1e-6

-- | Birth and death prior for bifurcating trees.
--
-- See Stadler, T., Mammalian phylogeny reveals recent diversification rate
-- shifts, Proceedings of the National Academy of Sciences, 108(15), 6187–6192
-- (2011). http://dx.doi.org/10.1073/pnas.1016876108.
--
-- The prior conditions on the time of origin.
--
-- TODO: Condition on time of origin and on survival.
--
-- TODO: Condition on time of origin and on the number of taxa.
--
-- XXX: This prior does not condition on survival.
--
-- XXX: This prior conditions on the origin, and not on the root node. This
-- affects: (1) the time, and (2) the split at the root, because the prior has
-- an additional multiplicative factor.
--
-- XXX: The prior does calculate the multiplicative combinatorial factor
-- relating the number of oriented labeled trees to the number of labeled trees
-- without orientation.
--
-- The aforementioned points are just multiplicative factors that don't
-- influence the stationary distribution of the MCMC run. However, they do
-- change the absolute value of the prior function, and hence, the posterior.
--
-- Call 'error' if
-- - The birth or death rate are negative.
-- - The sampling rate is zero or negative, or above 1.0.
-- - The tree is not bifurcating.
birthDeath ::
  -- | Birth rate.
  Double ->
  -- | Death rate.
  Double ->
  -- | Sampling rate.
  Double ->
  Tree Length a ->
  Log Double
birthDeath la mu rho
  | la < 0.0 = error "birthDeath: Birth rate is negative."
  | mu < 0.0 = error "birthDeath: Death rate is negative."
  | rho <= 0.0 = error "birthDeath: Sampling rate is zero or negative."
  | rho > 1.0 = error "birthDeath: Sampling rate is larger than 1.0."
  | epsNearCritical > abs (la - mu) = fst . birthDeathWith computeDENearCritical la mu rho (Exp $ log la)
  | otherwise = fst . birthDeathWith computeDE la mu rho (Exp $ log la)

birthDeathWith ::
  -- Computation of D and E. Set to normal or near critical formula.
  (Double -> Double -> Double -> Double -> (Double, Double)) ->
  -- Birth rate.
  Double ->
  -- Death rate.
  Double ->
  -- Sampling rate.
  Double ->
  -- Birth rate in log domain.
  Log Double ->
  Tree Length a ->
  (Log Double, Double)
birthDeathWith f la mu rho _ (Node br _ []) = first (Exp . log) $ f la mu (fromLength br) (1.0 - rho)
birthDeathWith f la mu rho logLa (Node br _ [l, r]) = (Exp (log dT) * dL * dR * logLa, eT)
  where
    (dL, eL) = birthDeathWith f la mu rho logLa l
    (dR, eR) = birthDeathWith f la mu rho logLa r
    (dT, eT) = f la mu (fromLength br) (eL * eR)
birthDeathWith _ _ _ _ _ _ = error "birthDeath: Tree is not bifurcating."

-- * Tests

--
-- >>> testTree1 :: Tree Double ()
-- >>> testTree1 = Node 1.0 () []
--
-- >>> birthDeath 1.2 3.2 testTree1
-- 5.8669248906043234e-2
--
-- >>> testTree2 :: Tree Double ()
-- >>> testTree2 = Node 0.0 () [Node 0.4 () [], Node 0.2 () [Node 0.2 () [], Node 0.2 () []]]
--
-- >>> birthDeath 1.2 3.2 testTree2
-- 3.978845396350806e-2

-- * Point process

--
-- There are differences in the conditions. The point process conditions on the
-- time of origin, and on the number of leaves. The dynamic programming approach
-- used above only conditions on the time of origin.
--
-- Question: Is the birth death distribution of the point process adequate, if
-- the topology is given?
--
-- How do we handle the critical, or nearly critical process? Is it fast enough
-- to combine the three distributions into one, and check during calculation of
-- the density? Or is it better to check for the three different possibilities
-- here, and then use the adequate distribution (status quo)?

-- -- Assume that node labels denote node heights.
-- birthDeathPointProcess ::
--   HasHeight a =>
--   -- | Birth rate.
--   Double ->
--   -- | Death rate.
--   Double ->
--   Tree Double a ->
--   Log Double
-- birthDeathPointProcess l m t
--   | l < 0.0 = error "birthDeath: Birth rate lambda is negative."
--   | m < 0.0 = error "birthDeath: Death rate mu is negative."
--   | otherwise = Exp $ logDensity d 0
--   where
--     d = BDD (label t + branch t) l m
