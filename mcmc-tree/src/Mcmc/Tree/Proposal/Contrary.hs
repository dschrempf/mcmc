-- |
-- Module      :  Mcmc.Tree.Proposal.Contrary
-- Description :  Contrary proposal between two trees
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Wed Mar  3 13:19:54 2021.
module Mcmc.Tree.Proposal.Contrary
  ( slideNodesAtContrarily,
    slideNodesContrarily,
    slideRootContrarily,
    scaleSubTreesAtContrarily,
    scaleSubTreesContrarily,
  )
where

import Control.Lens
import Control.Monad
import Data.Bifunctor
import ELynx.Tree
import Mcmc.Proposal
import Mcmc.Statistics.Types
import Mcmc.Tree.Lens
import Mcmc.Tree.Proposal.Common
import Mcmc.Tree.Proposal.Ultrametric
import Mcmc.Tree.Proposal.Unconstrained
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)

-- See also 'slideNodeAtUltrametricSimple'.
slideNodesAtContrarilySimple ::
  Path ->
  StandardDeviation Double ->
  TuningParameter ->
  ProposalSimple (HeightTree Double, Tree Double b)
slideNodesAtContrarilySimple pth sd t (tTr, rTr) g
  | null tTrChildren =
    error "slideNodesAtContrarilySimple: Sub tree of ultrametric tree is a leaf."
  | null rTrChildren =
    error "slideNodesAtContrarilySimple: Sub tree of unconstrained tree is a leaf."
  | otherwise = do
    (hTTrNode', q) <- truncatedNormalSample hTTrNode sd t hTTrOldestChild hTTrParent g
    -- Time tree.
    let setNodeHeight x = x & branchL .~ hTTrNode'
        tTr' = toTree $ modifyTree setNodeHeight tTrPos
    -- Rate tree.
    let -- Scaling factor of rate tree stem.
        xiStemR = if null pth then 1.0 else (hTTrParent - hTTrNode) / (hTTrParent - hTTrNode')
        -- Scaling factors of rate tree daughter branches excluding the stem.
        getXiR h = (hTTrNode - h) / (hTTrNode' - h)
        xisR = map getXiR hsTTrChildren
        scaleDaughterBranches (Node br lb trs) = Node br lb $ zipWith scaleUnconstrainedStem xisR trs
        -- If the root node is handled, do not scale the stem because no upper
        -- bound is set.
        f =
          if null pth
            then scaleDaughterBranches
            else scaleUnconstrainedStem xiStemR . scaleDaughterBranches
        rTr' = toTree $ modifyTree f rTrPos
    -- New state.
    let x' = (HeightTree tTr', rTr')
        -- jacobianTimeTree = Exp $ fromIntegral (nNodes - 1) * log xi
        -- jacobianRateTree = Exp $ fromIntegral (nBranches -1) * log xi' + log xiStem
        jacobian = Exp $ sum (map log xisR) + log xiStemR
    let
    return (x', q, jacobian)
  where
    -- Time tree.
    tTrPos = goPathUnsafe pth $ fromTree $ fromHeightTree tTr
    tTrFocus = current tTrPos
    tTrParent = current $ goParentUnsafe tTrPos
    hTTrNode = branch tTrFocus
    -- If the root node is handled, set the upper bound to +Infinity because no
    -- parent node exists.
    hTTrParent = if null pth then 1 / 0 else branch tTrParent
    tTrChildren = forest tTrFocus
    hsTTrChildren = map branch tTrChildren
    hTTrOldestChild = maximum hsTTrChildren
    -- Rate tree.
    rTrPos = goPathUnsafe pth $ fromTree rTr
    rTrFocus = current rTrPos
    rTrChildren = forest rTrFocus

-- | Slide nodes contrarily at given path.
--
-- This proposal acts on a pair of trees. The first tree is an ultrametric tree
-- (see "Mcmc.Tree.Proposal.Ultrametric"), usually the time tree with branches
-- denoting absolute or relative time. The second tree is an unconstrained tree
-- (see "Mcmc.Tree.Proposal.Unconstrained"), usually the rate tree with branches
-- denoting absolute or relative rates.
--
-- The specified nodes of both trees are slid contrarily. For example, if the
-- parent branch of the ultrametric tree node is shortened, the parent branch of
-- the unconstrained tree node is elongated. Ultrametricity is maintained.
--
-- A normal distribution truncated at the height of the parent node of the
-- ultrametric tree and the oldest daughter node of the ultrametric tree is used
-- to determine the new node height of the ultrametric tree. The rates are
-- changed accordingly.
--
-- See 'slideNodeAtUltrametric'.
--
-- NOTE: When applying to the root node (1) the tree heights change contrarily,
-- (2) no upper bound is used because no parent node exists, and (3) the stem of
-- the unconstrained tree is not changed because it does not map to a valid
-- branch of the ultrametric tree.
--
-- NOTE: Application of this proposal to trees with different topologies will
-- lead to unexpected behavior and possibly to run time errors.
--
-- Call 'error' if:
--
-- - The path is invalid.
--
-- - The path leads to a leaf.
--
-- - A node height or branch length is zero.
slideNodesAtContrarily ::
  -- | The topology of the tree is used to check the path.
  Tree e a ->
  Path ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree Double, Tree Double c)
slideNodesAtContrarily tr pth sd
  | not $ isValidPath tr pth = error $ "slideNodesAtContrarily: Path is invalid: " <> show pth <> "."
  | isLeafPath tr pth = error $ "slideNodesAtContrarily: Path leads to a leaf: " <> show pth <> "."
  | otherwise =
    createProposal
      description
      (slideNodesAtContrarilySimple pth sd)
      -- 1 for ultrametric node.
      -- 0 or 1 for unconstrained stem.
      -- n for unconstrained daughters.
      (PDimension $ 1 + nStem + nDaughters)
  where
    description = PDescription $ "Slide nodes contrarily; sd: " <> show sd
    nStem = if null pth then 0 else 1
    nDaughters = length $ forest $ current $ goPathUnsafe pth $ fromTree tr

-- | Slide the nodes of two trees contrarily.
--
-- See 'slideNodesAtContrarily'.
--
-- The weights are assigned as described in 'scaleSubTreesUltrametric'.
--
-- Do not scale the leaves.
slideNodesContrarily ::
  Tree e a ->
  HandleNode ->
  StandardDeviation Double ->
  PName ->
  -- | Minimum weight.
  PWeight ->
  -- | Maximum weight.
  PWeight ->
  Tune ->
  [Proposal (HeightTree Double, Tree Double c)]
slideNodesContrarily tr hn s n wMin wMax t =
  [ slideNodesAtContrarily tr pth s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      let focus = tr ^. subTreeAtL pth,
      let currentDepth = depth focus,
      -- Subtract 2 because leaves have depth one and are not scaled.
      let w = pWeight $ minimum [fromPWeight wMin + currentDepth - 2, fromPWeight wMax],
      -- Do not scale the leaves.
      not $ null $ forest focus,
      -- Filter other nodes.
      hn pth
  ]
  where
    name lb = n <> PName (" node " ++ show lb)

slideRootContrarilyJacobian ::
  -- Number of inner nodes.
  Int ->
  -- Scaling factor u for heights.
  Double ->
  -- Scaling factors xi_j for rates. Usually those are just two factors, but
  -- let's keep it general.
  [Double] ->
  Jacobian
slideRootContrarilyJacobian n u xis =
  Exp $
    sum $
      -- Minus n: Scaling the time tree node heights contrarily.
      fromIntegral (- n) * log u :
      -- Scaling the rate branches.
      map log xis

slideRootSimple ::
  Int ->
  StandardDeviation Double ->
  TuningParameter ->
  ProposalSimple (Double, HeightTree Double, Tree Double b)
slideRootSimple n s t (ht, HeightTree tTr, rTr) g = do
  let tTrHeight = branch tTr
  when
    (tTrHeight /= 1.0)
    ( error $
        "slideRootSimple: Height of relative time tree is different from 1.0: "
          <> show tTrHeight
          <> "."
    )
  -- Do not use 'genericContinuous' because the forward operator and the
  -- Jacobian function need access to the time tree node heights.
  (ht', q) <- truncatedNormalSample ht s t htOldestChild (1 / 0) g
  -- Scaling factor of absolute time tree height. This is the reverse scaling
  -- factor of the time tree node heights.
  let u = ht' / ht
  -- Scaling factors for rates.
  let getXi h = (1 - h) / (u - h)
      xis = map getXi htsChildren
  -- Compute new state.
  let tTr' = tTr & forestL %~ map (first (/ u))
      rTr' = rTr & forestL %~ zipWith scaleUnconstrainedStem xis
      j = slideRootContrarilyJacobian n u xis
      x' = (ht', HeightTree tTr', rTr')
  return (x', q, j)
  where
    htsChildren = map branch $ forest tTr
    -- Absolute height of oldest child.
    htOldestChild = ht * maximum htsChildren

-- | Specific proposal sliding the absolute time height while leaving the
-- absolute heights of internal time tree nodes untouched.
--
-- The proposal works on a parameter triple \((H, t, r)\), where \(H\) is the
-- absolute height of the time tree, \(t\) is a relative time tree, and \(r\) is
-- an absolute or relative rate tree.
--
-- Use a truncated normal distribution with given standard deviation to propose
-- a new height \(H'\) such that \(H'\) is larger than the highest daughter node
-- of the root. Let \(H'=H*u\). Scale all node heights of \(t\) contrarily. That
-- is, let \(I\) be the index set traversing the nodes of the time tree
-- excluding the leaves and the root. For any \(i \in I\), the node height
-- \(t_i\) will become \(t_i'=t_i/u\). Further, propose new rates for the rate
-- tree branches \(r_j\) leading to the root node. In particular,
-- \(r_j'=r_j\frac{1-t_j}{u-t_j}\), where \(t_j\) is the height of the node
-- corresponding to branch \(j\). In this way, the expected number of
-- substitutions on all branches stays constant.
--
-- Call 'error' if:
--
-- - The height of the relative time tree is different from 1.0.
slideRootContrarily ::
  -- | The topology of the tree is used to precompute the number of inner nodes.
  Tree e a ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Double, HeightTree Double, Tree Double c)
slideRootContrarily tr s =
  createProposal
    description
    (slideRootSimple n s)
    -- 1: Slide absolute time height.
    -- n: Scale inner nodes of time tree.
    -- k: Scale the rate tree branches leading to the root.
    (PDimension $ 1 + n + k)
  where
    description = PDescription $ "Slide root contrarily; sd: " ++ show s
    n = nInnerNodes tr
    k = length $ forest tr

-- See also 'scaleSubTreeAtUltrametricSimple'.
scaleSubTreeAtContrarilySimple ::
  -- Number of inner nodes.
  Int ->
  -- Number of branches.
  Int ->
  Path ->
  StandardDeviation Double ->
  TuningParameter ->
  ProposalSimple (HeightTree Double, Tree Double b)
scaleSubTreeAtContrarilySimple nNodes nBranches pth sd t (HeightTree tTr, rTr) g
  | null tTrChildren =
    error "scaleSubTreeAtContrarilySimple: Sub tree of ultrametric tree is a leaf."
  | null rTrChildren =
    error "scaleSubTreeAtContrarilySimple: Sub tree of unconstrained tree is a leaf."
  | otherwise = do
    (hTTrNode', q) <- truncatedNormalSample hTTrNode sd t 0 hTTrParent g
    let -- Scaling factor of time tree nodes heights (xi, not x_i).
        xiT = hTTrNode' / hTTrNode
        tTr' = toTree $ modifyTree (scaleUltrametricTreeF hTTrNode' xiT) tTrPos
    -- Rate tree.
    let -- Scaling factor of rate tree branches excluding the stem.
        xiR = recip xiT
        -- Scaling factor of rate tree stem.
        xiStemR = if null pth then 1.0 else (hTTrParent - hTTrNode) / (hTTrParent - hTTrNode')
        -- If the root node is handled, do not scale the stem because no upper
        -- bound is set.
        f =
          if null pth
            then scaleUnconstrainedTreeWithoutStemF xiR
            else scaleUnconstrainedStem xiStemR . scaleUnconstrainedTreeWithoutStemF xiR
        rTr' = toTree $ modifyTree f rTrPos
    -- New state.
    let x' = (HeightTree tTr', rTr')
        -- jacobianTimeTree = Exp $ fromIntegral (nNodes - 1) * log xi
        -- jacobianRateTree = Exp $ fromIntegral (nBranches -1) * log xi' + log xiStem
        jacobian = Exp $ fromIntegral (nNodes - nBranches) * log xiT + log xiStemR
    let
    return (x', q, jacobian)
  where
    -- Time tree.
    tTrPos = goPathUnsafe pth $ fromTree tTr
    tTrFocus = current tTrPos
    tTrParent = current $ goParentUnsafe tTrPos
    tTrChildren = forest tTrFocus
    hTTrNode = branch tTrFocus
    -- If the root node is handled, set the upper bound to +Infinity because no
    -- parent node exists.
    hTTrParent = if null pth then 1 / 0 else branch tTrParent
    -- Rate tree.
    rTrPos = goPathUnsafe pth $ fromTree rTr
    rTrFocus = current rTrPos
    rTrChildren = forest rTrFocus

-- | Scale the sub trees at given path.
--
-- This proposal acts on a pair of trees. The first tree is an ultrametric tree
-- (see "Mcmc.Tree.Proposal.Ultrametric"), usually the time tree with branches
-- denoting absolute or relative time. The second tree is an unconstrained tree
-- (see "Mcmc.Tree.Proposal.Unconstrained"), usually the rate tree with branches
-- denoting absolute or relative rates.
--
-- The sub trees of both trees are scaled contrarily. For example, if the sub
-- tree of the ultrametric tree is magnified, the sub tree of the unconstrained
-- tree is shrunk. In order to maintain ultrametricity of the ultrametric tree,
-- the stem of the sub tree is shortened (see 'scaleSubTreeAtUltrametric').
-- Correspondingly, the stem of the unconstrained tree is elongated.
--
-- A normal distribution truncated at the height of the parent node of the
-- ultrametric tree and the leaves is used to determine the new height of the
-- sub tree of the ultrametric tree.
--
-- For reference, please see 'scaleSubTreeAtUltrametric', and 'scaleTree'.
--
-- NOTE: When applying to the root node (1) the tree heights change contrarily,
-- (2) no upper bound is used because no parent node exists, and (3) the stem of
-- the unconstrained tree is not changed because it does not map to a valid
-- branch of the ultrametric tree.
--
-- NOTE: Application of this proposal to trees with different topologies will
-- lead to unexpected behavior and possibly to run time errors.
--
-- Call 'error' if:
--
-- - The path is invalid.
--
-- - The path leads to a leaf.
--
-- - A node height or branch length is zero.
scaleSubTreesAtContrarily ::
  -- | The topology of the tree is used to precompute the number of inner nodes.
  Tree e a ->
  Path ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree Double, Tree Double c)
scaleSubTreesAtContrarily tr pth sd
  | not $ isValidPath tr pth = error $ "scaleSubTreesAtContrarily: Path is invalid: " <> show pth <> "."
  | isLeafPath tr pth = error $ "scaleSubTreesAtContrarily: Path leads to a leaf: " <> show pth <> "."
  | otherwise =
    createProposal
      description
      (scaleSubTreeAtContrarilySimple nNodes nBranches pth sd)
      (PDimension $ nNodes + nBranches)
  where
    description = PDescription $ "Scale sub trees contrarily; sd: " <> show sd
    subtree = current $ goPathUnsafe pth $ fromTree tr
    nNodes = nInnerNodes subtree
    nBranches = length subtree

-- | Scale the sub trees of two trees contrarily.
--
-- See 'scaleSubTreesAtContrarily'.
--
-- The weights are assigned as described in 'scaleSubTreesUltrametric'.
--
-- Do not scale the leaves.
scaleSubTreesContrarily ::
  Tree e a ->
  HandleNode ->
  StandardDeviation Double ->
  PName ->
  -- | Minimum weight.
  PWeight ->
  -- | Maximum weight.
  PWeight ->
  Tune ->
  [Proposal (HeightTree Double, Tree Double c)]
scaleSubTreesContrarily tr hn s n wMin wMax t =
  [ scaleSubTreesAtContrarily tr pth s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      let focus = tr ^. subTreeAtL pth,
      let currentDepth = depth focus,
      -- Subtract 2 because leaves have depth one and are not scaled.
      let w = pWeight $ minimum [fromPWeight wMin + currentDepth - 2, fromPWeight wMax],
      -- Do not scale the leaves.
      not $ null $ forest focus,
      -- Filter other nodes.
      hn pth
  ]
  where
    name lb = n <> PName (" node " ++ show lb)
