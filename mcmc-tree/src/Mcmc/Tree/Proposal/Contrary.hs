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
  ( scaleSubTreeAtContrarily,
    scaleSubTreesContrarily,
  )
where

import Control.Lens
import Data.Bifunctor
import ELynx.Tree
import Mcmc.Proposal
import Mcmc.Statistics.Types
import Mcmc.Tree.Lens
import Mcmc.Tree.Proposal.Common
import Mcmc.Tree.Proposal.Ultrametric
import Mcmc.Tree.Types
import Numeric.Log

-- Scale all branches but the stem.
scaleRateTreeF :: Double -> Tree Length a -> Tree Length a
scaleRateTreeF xi tr = tr & forestL %~ map (first (lengthUnsafeL *~ xi))

scaleStem :: Double -> Tree Length a -> Tree Length a
scaleStem xi tr = tr & branchL . lengthUnsafeL *~ xi

-- See also 'scaleSubTreeAtUltrametricSimple'.
scaleSubTreeAtContrarilySimple ::
  -- Number of inner nodes.
  Int ->
  -- Number of branches.
  Int ->
  Path ->
  StandardDeviation ->
  TuningParameter ->
  ProposalSimple (HeightTree a, Tree Length b)
scaleSubTreeAtContrarilySimple nNodes nBranches pth sd t (tTr, rTr) g
  | null tTrChildren =
    error "scaleSubTreeAtContrarilySimple: Sub tree of time tree is a leaf."
  | null rTrChildren =
    error "scaleSubTreeAtContrarilySimple: Sub tree of rate tree is a leaf."
  | otherwise = do
    (hTTrNode', q) <- truncatedNormalSample hTTrNode sd t 0 hTTrParent g
    let -- Scaling factor of time tree branches excluding the stem (xi, not x_i).
        xi = hTTrNode' / hTTrNode
        tTr' = toTree $ modifyTree (scaleTreeF hTTrNode' xi) tTrPos
    -- Rate tree.
    let -- Scaling factor of rate tree branches excluding the stem.
        xi' = recip xi
        -- Scaling factor of rate tree stem.
        xiStem = (hTTrParent - hTTrNode') / (hTTrParent - hTTrNode)
        rTr' = toTree $ modifyTree (scaleStem xiStem . scaleRateTreeF xi') rTrPos
    -- New state.
    let x' = (tTr', rTr')
        -- jacobianTimeTree = Exp $ fromIntegral (nNodes - 1) * log xi
        -- jacobianRateTree = Exp $ fromIntegral (nBranches -1) * log xi' + log xiStem
        jacobian = Exp $ fromIntegral (nNodes - nBranches) * log xi + log xiStem
    let
    return (x', q, jacobian)
  where
    -- Time tree.
    tTrPos = goPathUnsafe pth $ fromTree tTr
    tTrFocus = current tTrPos
    tTrParent = current $ goParentUnsafe tTrPos
    tTrChildren = forest tTrFocus
    hTTrNode = fromHeight $ nodeHeight $ label tTrFocus
    hTTrParent = fromHeight $ nodeHeight $ label tTrParent
    -- Rate tree.
    rTrPos = goPathUnsafe pth $ fromTree rTr
    rTrFocus = current rTrPos
    rTrChildren = forest rTrFocus

-- | Scale the sub trees at given path.
--
-- This proposal acts on a pair of trees. The first tree is an ultrametric tree
-- (see "Mcmc.Tree.Proposal.Ultrametric"), usually the time tree with branches
-- denoting relative time. The second tree is an unconstrained tree (see
-- "Mcmc.Tree.Proposal.Unconstrained"), usually the rate tree with branches
-- denoting relative rates.
--
-- The sub trees of both trees are scaled contrarily. For example, if the sub
-- tree of the time tree is magnified, the sub tree of the rate tree is shrunk.
-- In order to maintain ultrametricity of the time tree, the stem of the sub
-- tree is shortened (see 'scaleSubTreeAtUltrametric'). Correspondingly, the stem
-- of the rate tree is elongated.
--
-- A normal distribution truncated at the height of the parent node and the
-- leaves is used to determine the new height of the sub tree.
--
-- For reference, please see 'scaleSubTreeAtUltrametric', and
-- 'Mcmc.Tree.Proposal.Unconstrained.scaleTree'.
--
-- NOTE: The topologies of the trees have to be equal. Different topologies will
-- lead to unexpected behavior and possibly to run time errors.
--
-- Call 'error' if:
--
-- - The path is invalid.
--
-- - The path is empty and leads to the root.
--
-- - The path leads to a leaf.
scaleSubTreeAtContrarily ::
  -- | The topology of the tree is used to precompute the number of inner nodes.
  --
  -- TODO: Use ELynx.Topology here.
  Tree e a ->
  Path ->
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree b, Tree Length c)
scaleSubTreeAtContrarily tr pth sd
  | null pth = error "scaleSubTreeAtContrarily: Path is empty."
  | not $ isValidPath tr pth = error $ "scaleSubTreeAtContrarily: Path is invalid: " <> show pth <> "."
  | isLeafPath tr pth = error $ "scaleSubTreeAtContrarily: Path leads to a leaf: " <> show pth <> "."
  | otherwise =
    createProposal
      description
      (scaleSubTreeAtContrarilySimple nNodes nBranches pth sd)
      (PDimension $ nNodes + nBranches)
  where
    description = PDescription $ "Scale sub trees contrarily; sd: " <> show sd
    subtree = current $ goPathUnsafe pth $ fromTree tr
    -- XXX: The number of nodes and branches could be calculated in one go.
    nNodes = nInnerNodes subtree
    nBranches = length subtree

-- | Scale the sub trees of the given trees constrarily.
--
-- See 'scaleSubTreeAtContrarily'.
--
-- Do not scale the root nor the leaves.
scaleSubTreesContrarily ::
  Tree e a ->
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  [Proposal (HeightTree b, Tree Length c)]
scaleSubTreesContrarily tr s n w t =
  [ scaleSubTreeAtContrarily tr pth s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      let focus = tr ^. subTreeAtUnsafeL pth,
      -- Do not scale the root.
      not $ null pth,
      -- Do not scale the leaves.
      not $ null $ forest focus
  ]
  where
    name lb = n <> PName (" node " ++ show lb)
