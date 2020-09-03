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
module Mcmc.Tree.Proposal
  ( -- * Nodes
    slideNodeUltrametric,

    -- * Branches
    scaleBranch,
    slideBranch,

    -- * Trees
    scaleTree,
    scaleTreeUltrametric,
    scaleSubTree,
    scaleSubTreeUltrametric,
    scaleTreesContrarily,
  )
where

import Control.Lens
import Data.Bifunctor
import ELynx.Data.Tree
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Proposal.Scale
import Mcmc.Proposal.Slide
import Mcmc.Tree.Lens
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Gamma
import Statistics.Distribution.TruncatedNormal
import System.Random.MWC

-- Minimum branch length.
eps :: Double
eps = 1e-8

-- -- Be careful, this will loop forever if the parameters are not well chosen.
-- truncatedNormal ::
--   PrimMonad m =>
--   Double ->
--   Double ->
--   Double ->
--   Double ->
--   Gen (PrimState m) ->
--   m Double
-- truncatedNormal a b m s g = do
--   x' <- normal m s g
--   case x' of
--     x
--       | x < a -> truncatedNormal a b m s g
--       | x > b -> truncatedNormal a b m s g
--       | otherwise -> return x

modifyBranch :: (e -> e) -> Tree e a -> Tree e a
modifyBranch f (Node br lb ts) = Node (f br) lb ts

-- A very specific function that samples a delta value from the truncated normal
-- distribution with given bounds [a,b] and also computes the required factor of
-- the Metropolis-Hastings proposal ratio.
truncatedNormalSample ::
  Double -> Double -> Double -> Double -> GenIO -> IO (Double, Log Double)
truncatedNormalSample ds t a b g = do
  let s = t * ds * (b - a)
      d = truncatedNormalDistr 0 s a b
  dx <- genContinuous d g
  -- Compute Metropolis-Hastings factor.
  let a' = a - dx
      b' = b - dx
      s' = t * ds * (b' - a')
      d' = truncatedNormalDistr 0 s' a' b'
      qXY = Exp $ logDensity d dx
      qYX = Exp $ logDensity d' (- dx)
  return (dx, qYX / qXY)

slideRootUltrametricSample ::
  Double ->
  Double ->
  Tree Double Double ->
  GenIO ->
  IO (Tree Double Double, Log Double)
slideRootUltrametricSample _ _ (Node _ _ []) _ =
  error "slideRootUltrametricSample: Cannot slide leaf node."
slideRootUltrametricSample ds t (Node br lb ts) g = do
  let br' = minimum $ map branch ts
      a = negate $ br - eps
      b = br' - eps
  (dx, q) <- truncatedNormalSample ds t a b g
  let tr' = Node (br + dx) (lb - dx) (map (modifyBranch (subtract dx)) ts)
  return (tr', q)

slideRootUltrametricSimple :: Double -> Double -> ProposalSimple (Tree Double Double)
slideRootUltrametricSimple ds t = ProposalSimple $ slideRootUltrametricSample ds t

-- | Slide a node.
--
-- A normal distribution truncated at the parent node and the closest daughter
-- node is used. The mean of the normal distribution is 0, the standard
-- deviation is determined by a given value multiplied with the domain of the
-- truncated distribution.
--
-- The node to slide is specified by a path.
--
-- Assume the branch and node labels denote branch length and node height,
-- respecitvely.
--
-- Call __error__ if given path is invalid or leads to a leaf.
slideNodeUltrametric ::
  -- | Path to node on tree.
  Path ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation multiplier.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
slideNodeUltrametric pth n w ds t = nodeAt pth @~ createProposal n w (slideRootUltrametricSimple ds) t

-- | Scale a branch.
--
-- The gamma distribution with given shape is used. The scale is set such that
-- the mean of the distribution is 1.0.
--
-- The branch to scale is specified by a path to a node.
scaleBranch ::
  -- | Path to node on tree.
  Path ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
scaleBranch pth n w s t = (nodeAt pth . rootBranch) @~ scaleUnbiased n w s t

-- | Slide a branch.
--
-- Use a normal distribution with mean 0 and given standard deviation.
--
-- The branch to slide is specified by a path to a node.
slideBranch ::
  -- | Path to node on tree.
  Path ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
slideBranch pth n w s t = (nodeAt pth . rootBranch) @~ slideSymmetric n w s t

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
-- Assume node labels denote node height.
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

-- | Scale all branches of sub tree induced by a given node with a gamma
-- distributed kernel of given shape. The scale is set such that the mean is
-- 1.0.
scaleSubTree ::
  -- | Path.
  Path ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
scaleSubTree pth n w k t = nodeAt pth @~ createProposal n w (scaleTreeSimple k) t

scaleRootUltrametricSample ::
  Double ->
  Double ->
  Tree Double Double ->
  GenIO ->
  IO (Tree Double Double, Log Double)
scaleRootUltrametricSample _ _ (Node _ _ []) _ =
  error "scaleRootUltrametricSample: Cannot scale sub tree of leaf node."
scaleRootUltrametricSample ds t (Node br lb ts) g = do
  let h = lb
      a = negate $ br - eps
      b = h - eps
  (dx, q) <- truncatedNormalSample ds t a b g
  let h' = lb - dx
      xi = h' / h
      tr' = Node (br + dx) h' $ map (bimap (* xi) (* xi)) ts
  return (tr', q)

scaleRootUltrametricSimple :: Double -> Double -> ProposalSimple (Tree Double Double)
scaleRootUltrametricSimple ds t = ProposalSimple $ scaleRootUltrametricSample ds t

-- | Slide the node at given path and scale the branches of the induced sub tree
-- so that the tree stays ultrametric.
--
-- An additive normal distribution truncated at the parent node (or the origin)
-- and the leaves is used to slide the given node. The mean of the normal
-- distribution is 0, the standard deviation is determined by a given value
-- multiplied with the domain of the truncated distribution.
--
-- Assume node labels denote node height.
--
-- Call __error__ if given path is invalid or leads to a leaf.
scaleSubTreeUltrametric ::
  -- | Path to node.
  Path ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation multiplier.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
scaleSubTreeUltrametric pth n w ds t = nodeAt pth @~ createProposal n w (scaleRootUltrametricSimple ds) t

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
--
-- XXX: For the first tree, assume that node labels denote node height.
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
