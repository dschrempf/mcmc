-- |
-- Module      :  Mcmc.Tree.Proposal.Ultrametric
-- Description :  Proposals preserving tree height
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Nov  4 11:48:56 2020.
--
-- For reasons of computational efficiency the functions working with
-- ultrametric trees use node labels directly storing node heights.
--
-- There is a distinction between proposals on unconstrained trees and
-- ultrametric trees, because the Jacobian matrices differ. For ultrametric
-- trees, some branch lengths are constrained, whereas for unconstrained trees,
-- all branch lengths can vary freely.
--
-- For the calculation of the Jacobian matrices of ultrametric trees, the
-- heights of the inner nodes are used as parameters.
--
-- Moreover, usage of proposals from "Mcmc.Tree.Proposal.Unconstrained" on
-- ultrametric trees will produce bogus trees with mismatching branch lengths
-- and node heights.
--
-- For proposals on unconstrained trees, see "Mcmc.Tree.Proposal.Unconstrained".
module Mcmc.Tree.Proposal.Ultrametric
  ( -- * Proposals
    slideNodeAtUltrametric,
    slideNodesUltrametric,
    scaleSubTreeAtUltrametric,
    scaleSubTreesUltrametric,
    pulleyUltrametric,

    -- * Helper functions
    nInnerNodes,
    scaleUltrametricTreeF,
  )
where

import Control.Exception
import Control.Lens hiding (children)
import Data.Bifunctor
import ELynx.Tree
import Mcmc.Proposal
import Mcmc.Statistics.Types
import Mcmc.Tree.Import ()
import Mcmc.Tree.Lens
import Mcmc.Tree.Proposal.Common
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)
import System.Random.MWC

slideNodeAtUltrametricSimple ::
  Path ->
  StandardDeviation Double ->
  TuningParameter ->
  ProposalSimple (HeightTree Double)
slideNodeAtUltrametricSimple pth s t tr g
  | null children = error "slideNodeAtUltrametricSimple: Cannot slide leaf."
  | otherwise = do
    (hNode', q) <- truncatedNormalSample hNode s t hChild hParent g
    let setNodeHeight x = x & branchL .~ hNode'
    -- The absolute value of the determinant of the Jacobian is 1.0.
    return (HeightTree $ toTree $ modifyTree setNodeHeight trPos, q, 1.0)
  where
    trPos = goPathUnsafe pth $ fromTree $ fromHeightTree tr
    focus = current trPos
    children = forest focus
    hNode = branch focus
    hChild = maximum $ map branch children
    parent = current $ goParentUnsafe trPos
    -- Set the upper bound to +Infinity if no parent node exists.
    hParent = if null pth then 1 / 0 else branch parent

-- | Slide node (for ultrametric trees).
--
-- For ultrametric trees, we cannot exclusively scale branches such as with
-- 'Mcmc.Tree.Proposal.Unconstrained.scaleBranch', because this would change the
-- heights of all descendant nodes. Consequently, if the proposal was used on a
-- non-root node, it would break ultrametricity of the tree. Instead, we can
-- slide node heights.
--
-- A normal distribution truncated at the heights of the parent node and the
-- closest child node is used.
--
-- NOTE: When sliding the root node the tree height changes and no upper bound
-- is used because no parent node exists.
--
-- Call 'error' if:
--
-- - The path is invalid.
--
-- - The path leads to a leaf.
slideNodeAtUltrametric ::
  -- | The topology of the tree is used to check the path.
  Tree e a ->
  -- | A zipper with given 'Path' has to be used for this proposal, because we
  -- need access to the parent.
  Path ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree Double)
slideNodeAtUltrametric tr pth ds
  | not $ isValidPath tr pth = error $ "slideNodeAtUltrametric: Path is invalid: " <> show pth <> "."
  | isLeafPath tr pth = error $ "slideNodeAtUltrametric: Path leads to a leaf: " <> show pth <> "."
  | otherwise = createProposal description (slideNodeAtUltrametricSimple pth ds) (PDimension 1)
  where
    description = PDescription $ "Slide node ultrametric; sd: " ++ show ds

-- | Slide the nodes of a given tree.
--
-- See 'slideNodeAtUltrametric'.
--
-- Do not slide leaves.
slideNodesUltrametric ::
  Tree e a ->
  HandleNode ->
  StandardDeviation Double ->
  -- | Base name of proposals.
  PName ->
  PWeight ->
  Tune ->
  [Proposal (HeightTree Double)]
slideNodesUltrametric tr hn s n w t =
  [ slideNodeAtUltrametric tr pth s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      -- Do not slide the leaves.
      not (isLeafPath tr pth),
      -- Filter other nodes.
      hn pth
  ]
  where
    name lb = n <> PName (" node " ++ show lb)

scaleSubTreeAtUltrametricSimple ::
  -- Number of inner nodes.
  Int ->
  Path ->
  StandardDeviation Double ->
  TuningParameter ->
  ProposalSimple (HeightTree Double)
scaleSubTreeAtUltrametricSimple n pth sd t tr g
  | null children = error "scaleSubTreeAtUltrametricSimple: Sub tree is a leaf."
  | otherwise = do
    (hNode', q) <- truncatedNormalSample hNode sd t 0 hParent g
    -- Scaling factor (xi, not x_i).
    let xi = hNode' / hNode
        -- (-1) because the root height has an additive change.
        jacobian = Exp $ fromIntegral (n - 1) * log xi
    return (HeightTree $ toTree $ modifyTree (scaleUltrametricTreeF hNode' xi) trPos, q, jacobian)
  where
    trPos = goPathUnsafe pth $ fromTree $ fromHeightTree tr
    focus = current trPos
    children = forest focus
    hNode = branch focus
    parent = current $ goParentUnsafe trPos
    -- Set the upper bound to +Infinity if no parent node exists.
    hParent = if null pth then 1 / 0 else branch parent

-- | Scale the node heights of the sub tree at given path.
--
-- A normal distribution with mean at the current node height, and truncated at
-- the height of the parent node and the leaves is used to determine the new
-- height of the sub tree.
--
-- NOTE: When scaling the root node the tree height changes and no upper bound
-- is used because no parent node exists.
--
-- Call 'error' if:
--
-- - The path is invalid.
--
-- - The path leads to a leaf.
scaleSubTreeAtUltrametric ::
  -- | The topology of the tree is used to precompute the number of inner nodes.
  Tree e a ->
  -- | A zipper with given 'Path' has to be used for this proposal, because we need
  -- access to the parent.
  Path ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree Double)
scaleSubTreeAtUltrametric tr pth sd
  | not $ isValidPath tr pth = error $ "scaleSubTreeAtUltrametric: Path is invalid: " <> show pth <> "."
  | isLeafPath tr pth = error $ "scaleSubTreeAtUltrametric: Path leads to a leaf: " <> show pth <> "."
  | otherwise =
    createProposal
      description
      (scaleSubTreeAtUltrametricSimple n pth sd)
      (PDimension n)
  where
    description = PDescription $ "Scale sub tree ultrametrc; sd: " ++ show sd
    n = nInnerNodes $ current $ goPathUnsafe pth $ fromTree tr

-- | Scale the sub trees of a given tree.
--
-- The proposal weights are set linearly according to the 'depth' of the node on
-- the tree. Minimum and maximum weights have to be provided.
--
-- See 'scaleSubTreeAtUltrametric'.
--
-- Do not scale the leaves.
scaleSubTreesUltrametric ::
  Tree e a ->
  HandleNode ->
  StandardDeviation Double ->
  -- | Base name of proposals.
  PName ->
  -- | Minimum weight at the leaves. If the minimum weight is larger than the
  -- maximum weight, the maximum weight will be assigned to all proposals.
  PWeight ->
  -- | Maximum weight for nodes closer to the root.
  PWeight ->
  Tune ->
  [Proposal (HeightTree Double)]
scaleSubTreesUltrametric tr hn s n wMin wMax t =
  [ scaleSubTreeAtUltrametric tr pth s (name lb) w t
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

-- See 'pulleyTruncatedNormalSample'. However, we have to honor more constraints
-- in the ultrametric case.
pulleyUltrametricTruncatedNormalSample ::
  StandardDeviation Double ->
  TuningParameter ->
  HeightTree Double ->
  GenIO ->
  IO (Double, Log Double)
pulleyUltrametricTruncatedNormalSample s t (HeightTree (Node ht _ [l, r]))
  | brL <= 0 =
    error $
      "pulleyUltrametricTruncatedNormalSample: Left branch is zero or negative: " ++ show brL ++ "."
  | brR <= 0 =
    error $
      "pulleyUltrametricTruncatedNormalSample: Right branch is zero or negative: " ++ show brR ++ "."
  | otherwise = truncatedNormalSample 0 s t a b
  where
    -- The new branch lengths are not allowed to exceed the height of the node.
    --
    -- Left and right branch length.
    brL = ht - branch l
    brR = ht - branch r
    -- The constraints are larger than 0.
    constraintRightBoundary = ht - brL
    constraintLeftBoundary = ht - brR
    a = negate $ minimum [brL, constraintLeftBoundary]
    b = minimum [brR, constraintRightBoundary]
pulleyUltrametricTruncatedNormalSample _ _ _ =
  error "pulleyUltrametricTruncatedNormalSample: Node is not bifurcating."

pulleyUltrametricSimple ::
  -- Number of inner nodes of left tree.
  Int ->
  -- Number of inner nodes of right tree.
  Int ->
  StandardDeviation Double ->
  TuningParameter ->
  ProposalSimple (HeightTree Double)
pulleyUltrametricSimple nL nR s t tr@(HeightTree (Node br lb [l, r])) g = do
  (u, q) <- pulleyUltrametricTruncatedNormalSample s t tr g
  -- Left.
  let hL = branch l
      hL' = hL - u
      -- Scaling factor left. (hL - u)/hL = (1.0 - u/hL).
      xiL = hL' / hL
  -- Right.
  let hR = branch r
      hR' = hR + u
      -- Scaling factor right. (hR + u)/hR = (1.0 + u/hR).
      xiR = hR' / hR
  let tr' = Node br lb [scaleUltrametricTreeF hL' xiL l, scaleUltrametricTreeF hR' xiR r]
  -- The derivation of the Jacobian matrix is very lengthy. Similar to before,
  -- we parameterize the right and left trees into the heights of all other
  -- internal nodes. However, the left and right node heights are now treated in
  -- a different way. For reference, I took a picture, 20201030_122839_DRO.jpg.
  --
  -- (-1) because the root height has an additive change.
  let jacobianL = Exp $ fromIntegral (nL - 1) * log xiL
      jacobianR = Exp $ fromIntegral (nR - 1) * log xiR
  return (HeightTree tr', q, jacobianL * jacobianR)
pulleyUltrametricSimple _ _ _ _ _ _ = error "pulleyUltrametricSimple: Node is not bifurcating."

-- | Use a node as a pulley.
--
-- See 'Mcmc.Tree.Proposal.Unconstrained.pulley' but for ultrametric trees. The
-- sub trees are scaled so that the tree remains ultrametric.
--
-- Call 'error' if:
--
-- - The node is not bifurcating.
--
-- - Left sub tree is a leaf.
--
-- - Right sub tree is a leaf.
pulleyUltrametric ::
  -- | The topology of the tree is used to precompute the number of inner nodes.
  Tree e a ->
  StandardDeviation Double ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree Double)
pulleyUltrametric (Node _ _ [l, r]) d
  | null (forest l) = error "pulleyUltrametric: Left sub tree is a leaf."
  | null (forest r) = error "pulleyUltrametric: Right sub tree is a leaf."
  | otherwise = createProposal description (pulleyUltrametricSimple nL nR d) (PDimension $ nL + nR)
  where
    description = PDescription $ "Pulley ultrametric; sd: " ++ show d
    nL = nInnerNodes l
    nR = nInnerNodes r
pulleyUltrametric _ _ = error "pulleyUltrametric: Node is not bifurcating."

-- | Calculate the number of inner nodes.
nInnerNodes :: Tree e a -> Int
nInnerNodes (Node _ _ []) = 0
nInnerNodes tr = 1 + sum (map nInnerNodes $ forest tr)

-- | A very specific function scaling an ultrametric tree.
scaleUltrametricTreeF ::
  -- | New root node height.
  Double ->
  -- | Scaling factor for other nodes. The scaling factor for inner node heights
  -- is also given, since it is calculated anyways by the calling functions.
  Double ->
  Tree Double Name ->
  Tree Double Name
scaleUltrametricTreeF h xi (Node _ lb ts) =
  Node h' lb $ map (first (* xi')) ts
  where
    xi' = assert (xi > 0) xi
    h' = assert (h > 0) h
