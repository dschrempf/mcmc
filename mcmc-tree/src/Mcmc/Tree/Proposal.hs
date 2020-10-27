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
-- For reasons of computational efficiency, functions working with ultrametric
-- trees (@functionNameUltrametric@) __assume the node labels denote node
-- height__ and handle these values accordingly.
module Mcmc.Tree.Proposal
  ( -- * Slide branches
    slideBranch,
    slideNodeUltrametric,

    -- * Scale trees
    scaleTree,
    scaleTreeUltrametric,
    scaleSubTreeUltrametric,
    scaleTreesContrarily,
    pulley,
    pulleyUltrametric,
  )
where

import Control.Lens
import Control.Monad
import Data.Bifunctor
import ELynx.Tree hiding (description)
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
  -- | Standard deviation.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
slideBranch s n w t = rootBranch @~ slideSymmetric s n w t

-- -- Minimum branch length.
-- eps :: Double
-- eps = 1e-12

-- A very specific function that samples a delta value from the truncated normal
-- distribution with given bounds [a,b] and also computes the required factor of
-- the Metropolis-Hastings proposal ratio.
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
  dx <- genContinuous d g
  -- Compute Metropolis-Hastings factor.
  let a' = a - dx
      b' = b - dx
      d' = either error id $ truncatedNormalDistr 0 s' a' b'
      qXY = Exp $ logDensity d dx
      qYX = Exp $ logDensity d' (- dx)
  return (dx, qYX / qXY)

-- The branch is elongated by dx. So if dx is positive, the node height is
-- reduced.
slideNodeUltrametricF :: Double -> Tree Double Double -> Tree Double Double
slideNodeUltrametricF dx (Node br h ts) = Node (br + dx) (h - dx) (map (applyStem (subtract dx)) ts)

slideNodeUltrametricSimple ::
  -- Standard deviation.
  Double ->
  -- Tuning parameter.
  Double ->
  ProposalSimple (Tree Double Double)
slideNodeUltrametricSimple _ _ (Node _ _ []) _ =
  error "slideNodeUltrametricSample: Cannot slide leaf node."
slideNodeUltrametricSimple s t tr@(Node br _ ts) g = do
  when
    (br <= 0)
    (error $ "slideNodeUltrametricSimple: Parent branch length is zero or negative: " ++ show br ++ ".")
  when
    (br' <= 0)
    (error $ "slideNodeUltrametricSimple: Minimum branch length is zero or negative: " ++ show br' ++ ".")
  (dx, q) <- truncatedNormalSample s t a b g
  return (slideNodeUltrametricF dx tr, q)
  where
    br' = minimum $ map branch ts
    -- a = negate $ br - eps
    a = negate br
    -- b = br' - eps
    b = br'

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
  -- | Standard deviation.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
slideNodeUltrametric ds = createProposal description (slideNodeUltrametricSimple ds)
  where
    description = "Slide node ultrametric; sd: " ++ show ds

scaleTreeSimple :: Double -> Double -> ProposalSimple (Tree Double a)
scaleTreeSimple k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (\tr x -> first (* x) tr)
    (Just recip)

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
scaleTree ::
  -- | Shape.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
scaleTree k = createProposal description (scaleTreeSimple k)
  where
    description = "Scale tree; shape: " ++ show k

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
  -- | Shape.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
scaleTreeUltrametric k = createProposal description (scaleTreeUltrametricSimple k)
  where
    description = "Scale tree ultrametric; shape: " ++ show k

-- The branch is elongated by dx. So if dx is positive, the node height is
-- reduced.
slideBranchScaleSubTreeF :: Double -> Tree Double Double -> Tree Double Double
slideBranchScaleSubTreeF dx (Node br h ts) = Node (br + dx) h' $ map (bimap (* xi) (* xi)) ts
  where
    h' = h - dx
    xi = h' / h

scaleSubTreeUltrametricSimple ::
  Double ->
  Double ->
  ProposalSimple (Tree Double Double)
scaleSubTreeUltrametricSimple _ _ (Node _ _ []) _ =
  error "scaleSubTreeUltrametricSample: Cannot scale sub tree of leaf node."
scaleSubTreeUltrametricSimple ds t tr g = do
  when
    (br <= 0)
    (error $ "scaleSubTreeUltrametricSimple: Parent branch length is zero or negative: " ++ show br ++ ".")
  when
    (ht <= 0)
    (error $ "scaleSubTreeUltrametricSimple: Node height is zero or negative: " ++ show ht ++ ".")
  (dx, q) <- truncatedNormalSample ds t a b g
  return (slideBranchScaleSubTreeF dx tr, q)
  where
    br = branch tr
    ht = label tr
    -- a = negate $ branch tr - eps
    a = negate br
    -- b = label tr - eps
    b = ht

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
  -- | Standard deviation.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
scaleSubTreeUltrametric sd = createProposal description (scaleSubTreeUltrametricSimple sd)
  where
    description = "Scale subtree ultrametrc; sd: " ++ show sd

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
  -- | Shape.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double, Tree Double a)
scaleTreesContrarily k = createProposal description (scaleTreesContrarilySimple k)
  where
    description = "Scale trees contrarily; shape: " ++ show k

-- See 'truncatedNormalSample'. Dx is added to the left branch. I.e., if dx is
-- positive, the left branch is elongated.
pulleyTruncatedNormalSample :: Double -> Double -> Tree Double a -> GenIO -> IO (Double, Log Double)
pulleyTruncatedNormalSample s t (Node _ _ [l, r]) = do
  when
    (brL <= 0)
    (error $ "pulleyTruncatedNormalSample: Left branch is zero or negative: " ++ show brL ++ ".")
  when
    (brR <= 0)
    (error $ "pulleyTruncatedNormalSample: Right branch is zero or negative: " ++ show brR ++ ".")
  truncatedNormalSample s t a b
  where
    brL = branch l
    brR = branch r
    -- a = negate $ branch l - eps
    a = negate brL
    -- b = branch r - eps
    b = brR
pulleyTruncatedNormalSample _ _ _ = error "pulleyTruncatedNormalSample: Node is not bifurcating."

pulleySimple :: Double -> Double -> ProposalSimple (Tree Double a)
pulleySimple s t tr@(Node br lb [l, r]) g = do
  (dx, q) <- pulleyTruncatedNormalSample s t tr g
  let tr' = Node br lb [applyStem (+ dx) l, applyStem (subtract dx) r]
  return (tr', q)
pulleySimple _ _ _ _ = error "pulleySimple: Node is not bifurcating."

-- | Use a node as a pulley.
--
-- For bifurcating nodes; change the daughter branch lengths contrarily. The
-- first and second daughter branches are elongated/shortened by the same
-- amount.
--
-- The height of the two sub trees changes.
pulley ::
  -- | Standard deviation.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
pulley s = createProposal description (pulleySimple s)
  where
    description = "Pulley; sd: " ++ show s

pulleyUltrametricSimple :: Double -> Double -> ProposalSimple (Tree Double Double)
pulleyUltrametricSimple s t tr@(Node br lb [l, r]) g = do
  (dx, q) <- pulleyTruncatedNormalSample s t tr g
  let tr' = Node br lb [slideBranchScaleSubTreeF dx l, slideBranchScaleSubTreeF (negate dx) r]
  return (tr', q)
pulleyUltrametricSimple _ _ _ _ = error "pulleyUltrametricSimple: Node is not bifurcating."

-- | Use a node as a pulley.
--
-- See 'pulley', but for ultrametric trees. The sub trees are scaled such that
-- the tree heights are conserved and the tree remains ultrametric.
--
-- __Assume the node labels denote node height__.
pulleyUltrametric ::
  -- | Standard deviation.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double Double)
pulleyUltrametric d = createProposal description (pulleyUltrametricSimple d)
  where
    description = "Pulley ultrametric; sd: " ++ show d
