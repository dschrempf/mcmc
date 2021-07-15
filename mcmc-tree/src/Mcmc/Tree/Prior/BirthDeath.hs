-- |
-- Module      :  Mcmc.Tree.Prior.BirthDeath
-- Description :  Calculate probability of a tree assuming birth and death process
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Aug  4 20:37:11 2020.
module Mcmc.Tree.Prior.BirthDeath
  ( BirthRate,
    DeathRate,
    SamplingRate,
    ConditionOn (..),
    birthDeath,
  )
where

-- I tested these functions many times, and obtained correct results! I also
-- compared the results to RevBayes, see below.

import ELynx.Tree
import Mcmc.Prior
import Numeric.Log

-- | Type synonym to indicate a birth rate.
type BirthRate a = a

-- | Type synonym to indicate a death rate.
type DeathRate a = a

-- | Type synonym to indicate a sampling rate.
type SamplingRate a = a

-- Compute probabilities D and E at the top of the branch.
--
-- Stadler, T. Mammalian phylogeny reveals recent diversification rate shifts.
-- Proceedings of the National Academy of Sciences 108, 6187–6192 (2011).
--
-- E is given in Eq. [1].
--
-- D is given in Eq. [2] without the multiplicative value coming from the
-- boundary conditions.
--
-- Correct results:
-- >>> computeDE 1.2 3.2 1.0 0.3
-- (7.283127121752474e-2,0.9305035687810801)
computeDE ::
  RealFloat a =>
  -- Lambda.
  BirthRate a ->
  -- Mu.
  DeathRate a ->
  -- Rho.
  SamplingRate a ->
  -- Branch length.
  a ->
  -- E at bottom of branch. Storage in log domain not necessary.
  a ->
  -- D, E. Storage in log domain not necessary.
  (a, a)
computeDE la mu rho dt e0 = (nomD / denom / denom, nomE / denom)
  where
    d = la - mu
    x = exp (- d * dt)
    c = (1.0 - rho) + rho * e0
    y = (mu - c * la) * x
    nomD = d * d * x
    c' = c - 1.0
    nomE = mu * c' + y
    denom = la * c' + y
{-# INLINE computeDE #-}

-- Compute probabilities D and E at the top of the branch for la ~= mu.
computeDENearCritical ::
  RealFloat a =>
  -- Lambda.
  BirthRate a ->
  -- Mu.
  DeathRate a ->
  -- Rho.
  SamplingRate a ->
  -- Branch length.
  a ->
  -- E at bottom of branch. Storage in log domain not necessary.
  a ->
  -- D, E. Storage in log domain not necessary.
  (a, a)
computeDENearCritical la mu rho dt e0 = (nomD / denom / denom, nomE / denom)
  where
    d = la - mu
    c = (1.0 - rho) + rho * e0
    y = (mu - c * la) * dt
    nomD = 1 - d * dt
    nomE = c + y
    denom = 1.0 + y
{-# INLINE computeDENearCritical #-}

-- Require near critical process if birth and death rates are closer than this value.
epsNearCritical :: Fractional a => a
epsNearCritical = 1e-6

-- | Condition on the time of origin or the time of the most recent common
-- ancestor (MRCA).
--
-- If 'ConditionOnTimeOfOrigin' is used, the stem needs to be strictly positive.
--
-- NOTE: It is not yet possible to condition on:
--
-- - the time of origin and survival (or MRCA and survival);
--
-- - the time of origin and the number of taxa (or MRCA and the number of taxa).
data ConditionOn = ConditionOnTimeOfOrigin | ConditionOnTimeOfMrca

-- | Birth and death prior for bifurcating trees.
--
-- See Stadler, T., Mammalian phylogeny reveals recent diversification rate
-- shifts, Proceedings of the National Academy of Sciences, 108(15), 6187–6192
-- (2011). http://dx.doi.org/10.1073/pnas.1016876108.
--
-- NOTE: The prior __does not calculate the multiplicative combinatorial
-- factor__ relating the number of oriented labeled trees to the number of
-- labeled trees without orientation. This is an issue if the prior is
-- calculated for trees with different topologies!
--
-- Call 'error' if
--
-- - The birth or death rate are negative.
--
-- - The sampling rate is zero or negative, or above 1.0.
--
-- - The tree is not bifurcating.
birthDeath ::
  RealFloat a =>
  ConditionOn ->
  BirthRate a ->
  DeathRate a ->
  SamplingRate a ->
  PriorFunctionG (Tree a b) a
birthDeath ConditionOnTimeOfOrigin la mu rho t
  | la < 0.0 = error "birthDeath: Birth rate is negative."
  | mu < 0.0 = error "birthDeath: Death rate is negative."
  | rho <= 0.0 = error "birthDeath: Sampling rate is zero or negative."
  | rho > 1.0 = error "birthDeath: Sampling rate is larger than 1.0."
  | epsNearCritical > abs (la - mu) = fst $ birthDeathWith computeDENearCritical la mu rho t
  | otherwise = fst $ birthDeathWith computeDE la mu rho t
birthDeath ConditionOnTimeOfMrca la mu rho (Node _ _ [l, r]) =
  birthDeath ConditionOnTimeOfOrigin la mu rho l * birthDeath ConditionOnTimeOfOrigin la mu rho r
birthDeath ConditionOnTimeOfMrca _ _ _ _ =
  error "birthDeath: Tree is not bifurcating."
{-# SPECIALIZE birthDeath ::
  ConditionOn ->
  Double ->
  Double ->
  Double ->
  PriorFunction (Tree Double b)
  #-}

birthDeathWith ::
  RealFloat a =>
  -- Computation of D and E. Set to normal or near critical formula.
  (BirthRate a -> DeathRate a -> SamplingRate a -> a -> a -> (a, a)) ->
  BirthRate a ->
  DeathRate a ->
  SamplingRate a ->
  Tree a b ->
  -- Return (log D, E).
  (Log a, a)
-- First case of the boundary conditions given after Eq. [4].
birthDeathWith f la mu rho (Node br _ [l, r]) = (Exp (log (dT * la)) * dL + dR, eT)
  where
    (dL, eL) = birthDeathWith f la mu rho l
    -- (dR, eR) = birthDeathWith f la mu rho r
    --
    -- XXX: eL and eR should be the same. However, there will be hardly any
    -- speed benefit because D and E are calculated simultaneously. Should I
    -- separate the calculation of D and E?
    --
    -- XXX: I should actually check that eL and eR indeed are the same, but that
    -- would also slow down the calculation.
    (dR, _) = birthDeathWith f la mu rho r
    -- D and E at the top of the internal branch. Since we are treating internal
    -- nodes here, we use rho=1.0. In the future, one may also allow values
    -- below 1.0 modeling, for example, catastrophes killing a certain
    -- percentage of all living species.
    (dT, eT) = f la mu 1.0 br eL
-- Second case of the boundary conditions given after Eq. [4].
birthDeathWith f la mu rho (Node br _ [c]) = (Exp (log (dT * rho)) * d, eT)
  where
    (d, e) = birthDeathWith f la mu rho c
    (dT, eT) = f la mu 1.0 br e
-- Third case of the boundary conditions given after Eq. [4].
birthDeathWith f la mu rho (Node br _ []) = (Exp $ log $ dT * rho, eT)
  where
    -- D and E at the top of the external branch. We use the given sampling
    -- probability here.
    (dT, eT) = f la mu rho br 0
birthDeathWith _ _ _ _ _ = error "birthDeathWith: Tree is multifurcating."

-- * Tests

--
-- >>> let testTree1 = Node 1.0 () [] :: Tree Length ()
-- >>> birthDeath WithStem 1.2 3.2 1.0 testTree1
-- 5.8669248906043234e-2
--
-- >>> let testTree2 = Node 0.0 () [Node 0.4 () [], Node 0.2 () [Node 0.2 () [], Node 0.2 () []]] :: Tree Length ()
-- >>> birthDeath WithStem 1.2 3.2 1.0 testTree2
-- 4.3357752474276125e-2
--
-- The following computations are checked against RevBayes:
--
-- >>> let t = parseNewick Standard "(((a:1.0,b:1.0):1.0,c:2.0):1.0,d:3.0):0.0;"
--
-- >>> map (\mu -> ln $ 1/3* (birthDeath WithoutStem 1.0 mu 1.0 $ either error id $ phyloToLengthTree t)) [0, 0.01, 0.05, 0.1, 0.2, 0.5]
-- [-10.09861228866811,-10.07675364864067,-9.993307032921498,-9.898174270006024,-9.73975910235509,-9.54137886890279]
--
-- >>> map (\rho -> ln $ 1/3* (birthDeath WithoutStem 1.0 0.0 rho $ either error id $ phyloToLengthTree t)) [1.0, 0.9, 0.8]
-- [-10.09861228866811,-9.809211822253452,-9.498032504556043]
--
-- >>> ln $ 1/3 * (birthDeath WithoutStem 0.2 0.5 0.8 $ either error id $ phyloToLengthTree t)
-- -9.700151607658995

-- * Point process

--
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
--   HasLength a =>
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
