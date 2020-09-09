-- TODO: Proposals on tree topologies.
-- - NNI
-- - Narrow (what is this, see RevBayes)
-- - FNPR (dito)

-- |
-- Module      :  Mcmc.Tree.Proposal
-- Description :  Proposals on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 09:10:07 2020.
--
-- For reasons of computational efficiency, all functions named
-- @nameUltrametric@ __assume the node labels denote node height__ and handle
-- these values accordingly.
module Mcmc.Tree.Proposal
  ( -- * Slide branches
    slideBranch,
    slideNodeUltrametric,

    -- * Scale trees
    scaleTree,
    scaleTreeUltrametric,
    scaleSubTreeUltrametric,
    scaleTreesContrarily,
  )
where

import Control.Lens
import Data.Bifunctor
import ELynx.Tree
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Proposal.Slide
import Mcmc.Tree.Lens
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Gamma
import Statistics.Distribution.TruncatedNormal
import System.Random.MWC

-- | Slide branch.
--
-- Use a normal distribution with mean 0 and given standard deviation.
--
-- This proposal slides the root branch. To slide other branches, see 'nodeAt'.
-- For example, @nodeAt path @~ slideNodeUltrametric ...@.
slideBranch ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
slideBranch n w s t = rootBranch @~ slideSymmetric n w s t

-- Minimum branch length.
eps :: Double
eps = 1e-8

-- A very specific function that samples a delta value from the truncated normal
-- distribution with given bounds [a,b] and also computes the required factor of
-- the Metropolis-Hastings proposal ratio.
truncatedNormalSample ::
  Double -> Double -> Double -> Double -> GenIO -> IO (Double, Log Double)
truncatedNormalSample s t a b g = do
  let s' = t * s
      d = truncatedNormalDistr 0 s' a b
  dx <- genContinuous d g
  -- Compute Metropolis-Hastings factor.
  let a' = a - dx
      b' = b - dx
      d' = truncatedNormalDistr 0 s' a' b'
      qXY = Exp $ logDensity d dx
      qYX = Exp $ logDensity d' (- dx)
  return (dx, qYX / qXY)

slideNodeUltrametricSample ::
  Double ->
  Double ->
  Tree Double Double ->
  GenIO ->
  IO (Tree Double Double, Log Double)
slideNodeUltrametricSample _ _ (Node _ _ []) _ =
  error "slideNodeUltrametricSample: Cannot slide leaf node."
slideNodeUltrametricSample ds t (Node br lb ts) g = do
  let br' = minimum $ map branch ts
      a = negate $ br - eps
      b = br' - eps
  (dx, q) <- truncatedNormalSample ds t a b g
  let tr' = Node (br + dx) (lb - dx) (map (applyStem (subtract dx)) ts)
  return (tr', q)

slideNodeUltrametricSimple :: Double -> Double -> ProposalSimple (Tree Double Double)
slideNodeUltrametricSimple ds t = ProposalSimple $ slideNodeUltrametricSample ds t

-- | Slide node (for ultrametric trees).
--
-- For ultrametric trees, we cannot exclusively slide the branch such as with
-- 'slideBranch', because this would change the height and if the proposal is
-- used on a non-root node, it would break ultrametricity of the tree. Instead,
-- we need to slide the root node. That is, when the stem is elongated, we need
-- to shorten the daughter branches, and vice versa, such that the tree height
-- is conserved.
--
-- A normal distribution truncated at the origin and the closest daughter node
-- is used.
--
-- __Assume the node labels denote node height__.
slideNodeUltrametric ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
slideNodeUltrametric n w ds = createProposal n w (slideNodeUltrametricSimple ds)

scaleTreeSimple :: Double -> Double -> ProposalSimple (Tree Double a)
scaleTreeSimple k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (\tr x -> first (* x) tr)
    (Just recip)

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
scaleTree ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
scaleTree n w k = createProposal n w (scaleTreeSimple k)

scaleTreeUltrametricSimple :: Double -> Double -> ProposalSimple (Tree Double Double)
scaleTreeUltrametricSimple k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (\tr x -> bimap (* x) (* x) tr)
    (Just recip)

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
--
-- __Assume node labels denote node height__.
--
-- __The height is changed.__ Do not use this proposal on a sub tree of an
-- ultrametric tree. Instead, use 'scaleSubTreeUltrametric'.
scaleTreeUltrametric ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
scaleTreeUltrametric n w k = createProposal n w (scaleTreeUltrametricSimple k)

scaleSubTreeUltrametricSample ::
  Double ->
  Double ->
  Tree Double Double ->
  GenIO ->
  IO (Tree Double Double, Log Double)
scaleSubTreeUltrametricSample _ _ (Node _ _ []) _ =
  error "scaleSubTreeUltrametricSample: Cannot scale sub tree of leaf node."
scaleSubTreeUltrametricSample ds t (Node br lb ts) g = do
  let h = lb
      a = negate $ br - eps
      b = h - eps
  (dx, q) <- truncatedNormalSample ds t a b g
  let h' = lb - dx
      xi = h' / h
      tr' = Node (br + dx) h' $ map (bimap (* xi) (* xi)) ts
  return (tr', q)

scaleSubTreeUltrametricSimple :: Double -> Double -> ProposalSimple (Tree Double Double)
scaleSubTreeUltrametricSimple ds t = ProposalSimple $ scaleSubTreeUltrametricSample ds t

-- | Scale the branches of the sub tree and slide the root branch so that the
-- tree height is conserved.
--
-- See also 'scaleTreeUltrametric'.
--
-- A normal distribution truncated at the parent node (or the origin) and the
-- leaves is used to slide the given node.
--
-- __Assume the node labels denote node height__.
scaleSubTreeUltrametric ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
scaleSubTreeUltrametric n w ds = createProposal n w (scaleSubTreeUltrametricSimple ds)

contra :: (Tree Double Double, Tree Double a) -> Double -> (Tree Double Double, Tree Double a)
contra (s, t) x = (bimap (* x) (* x) s, first (/ x) t)

scaleTreesContrarilySimple :: Double -> Double -> ProposalSimple (Tree Double Double, Tree Double a)
scaleTreesContrarilySimple k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    contra
    (Just recip)

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
--
-- The two trees are scaled contrarily so that the product of their heights
-- stays constant. Contrary proposals are useful when parameters are confounded.
scaleTreesContrarily ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double, Tree Double a)
scaleTreesContrarily n w k = createProposal n w (scaleTreesContrarilySimple k)
