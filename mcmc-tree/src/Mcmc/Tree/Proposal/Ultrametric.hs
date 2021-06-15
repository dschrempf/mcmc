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
  StandardDeviation ->
  TuningParameter ->
  ProposalSimple (HeightTree a)
slideNodeAtUltrametricSimple pth s t tr g
  | null children = error "slideNodeAtUltrametricSimple: Cannot slide leaf."
  | otherwise = do
    (hNode', q) <- truncatedNormalSample hNode s t hChild hParent g
    let setNodeHeight x =
          x & labelL . nodeHeightL . heightL .~ hNode'
    -- The absolute value of the determinant of the Jacobian is 1.0.
    return (toTree $ modifyTree setNodeHeight trPos, q, 1.0)
  where
    trPos = goPathUnsafe pth $ fromTree tr
    focus = current trPos
    children = forest focus
    hNode = fromHeight $ nodeHeight $ label focus
    hChild = fromHeight $ maximum $ map (nodeHeight . label) children
    parent = current $ goParentUnsafe trPos
    -- Set the upper bound to +Infinity if no parent node exists.
    hParent = if null pth then 1 / 0 else fromHeight $ nodeHeight $ label parent

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
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree b)
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
  HandleLayer ->
  StandardDeviation ->
  -- | Base name of proposals.
  PName ->
  PWeight ->
  Tune ->
  [Proposal (HeightTree b)]
slideNodesUltrametric tr hd s n w t =
  [ slideNodeAtUltrametric tr pth s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      -- Do not slide the leaves.
      not (isLeafPath tr pth),
      -- Filter other layers.
      hd $ length pth
  ]
  where
    name lb = n <> PName (" node " ++ show lb)

scaleSubTreeAtUltrametricSimple ::
  -- Number of inner nodes.
  Int ->
  Path ->
  StandardDeviation ->
  TuningParameter ->
  ProposalSimple (HeightTree a)
scaleSubTreeAtUltrametricSimple n pth sd t tr g
  | null children = error "scaleSubTreeAtUltrametricSimple: Sub tree is a leaf."
  | otherwise = do
    (hNode', q) <- truncatedNormalSample hNode sd t 0 hParent g
    -- Scaling factor (xi, not x_i).
    let xi = hNode' / hNode
        -- (-1) because the root height has an additive change.
        jacobian = Exp $ fromIntegral (n - 1) * log xi
    return (toTree $ modifyTree (scaleUltrametricTreeF hNode' xi) trPos, q, jacobian)
  where
    trPos = goPathUnsafe pth $ fromTree tr
    focus = current trPos
    children = forest focus
    hNode = fromHeight $ nodeHeight $ label focus
    parent = current $ goParentUnsafe trPos
    -- Set the upper bound to +Infinity if no parent node exists.
    hParent = if null pth then 1 / 0 else fromHeight $ nodeHeight $ label parent

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
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree b)
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
-- NOTE: The layer (see 'HandleLayer') and the 'depth' of a node are somewhat
-- opposing concepts.
--
-- See 'scaleSubTreeAtUltrametric'.
--
-- Do not scale the leaves.
scaleSubTreesUltrametric ::
  Tree e a ->
  HandleLayer ->
  StandardDeviation ->
  -- | Base name of proposals.
  PName ->
  -- | Minimum weight at the leaves. If the minimum weight is larger than the
  -- maximum weight, the maximum weight will be assigned to all proposals.
  PWeight ->
  -- | Maximum weight for nodes closer to the root.
  PWeight ->
  Tune ->
  [Proposal (HeightTree b)]
scaleSubTreesUltrametric tr hd s n wMin wMax t =
  [ scaleSubTreeAtUltrametric tr pth s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      let focus = tr ^. subTreeAtL pth,
      let currentLayer = length pth,
      let currentDepth = depth focus,
      -- Subtract 2 because leaves have depth one and are not scaled.
      let w = pWeight $ minimum [fromPWeight wMin + currentDepth - 2, fromPWeight wMax],
      -- Do not scale the leaves.
      not $ null $ forest focus,
      -- Filter other layers.
      hd currentLayer
  ]
  where
    name lb = n <> PName (" node " ++ show lb)

-- See 'pulleyTruncatedNormalSample'. However, we have to honor more constraints
-- in the ultrametric case.
pulleyUltrametricTruncatedNormalSample ::
  StandardDeviation ->
  TuningParameter ->
  HeightTree a ->
  GenIO ->
  IO (Double, Log Double)
pulleyUltrametricTruncatedNormalSample s t (Node _ lb [l, r])
  | brL <= 0 =
    error $
      "pulleyUltrametricTruncatedNormalSample: Left branch is zero or negative: " ++ show brL ++ "."
  | brR <= 0 =
    error $
      "pulleyUltrametricTruncatedNormalSample: Right branch is zero or negative: " ++ show brR ++ "."
  | otherwise = truncatedNormalSample 0 s t a b
  where
    -- The new branch lengths are not allowed to exceed the height of the node.
    ht = nodeHeight lb
    -- Left and right branch length.
    brL = ht - nodeHeight (label l)
    brR = ht - nodeHeight (label r)
    -- The constraints are larger than 0.
    constraintRightBoundary = ht - brL
    constraintLeftBoundary = ht - brR
    a = fromHeight $ negate $ minimum [brL, constraintLeftBoundary]
    b = fromHeight $ minimum [brR, constraintRightBoundary]
pulleyUltrametricTruncatedNormalSample _ _ _ =
  error "pulleyUltrametricTruncatedNormalSample: Node is not bifurcating."

pulleyUltrametricSimple ::
  -- Number of inner nodes of left tree.
  Int ->
  -- Number of inner nodes of right tree.
  Int ->
  StandardDeviation ->
  TuningParameter ->
  ProposalSimple (HeightTree a)
pulleyUltrametricSimple nL nR s t tr@(Node br lb [l, r]) g = do
  (u, q) <- pulleyUltrametricTruncatedNormalSample s t tr g
  -- Left.
  let hL = nodeHeight $ label l
      hL' = fromHeight hL - u
      -- Scaling factor left. (hL - u)/hL = (1.0 - u/hL).
      xiL = hL' / fromHeight hL
  -- Right.
  let hR = nodeHeight $ label r
      hR' = fromHeight hR + u
      -- Scaling factor right. (hR + u)/hR = (1.0 + u/hR).
      xiR = hR' / fromHeight hR
  let tr' = Node br lb [scaleUltrametricTreeF hL' xiL l, scaleUltrametricTreeF hR' xiR r]
  -- The derivation of the Jacobian matrix is very lengthy. Similar to before,
  -- we parameterize the right and left trees into the heights of all other
  -- internal nodes. However, the left and right node heights are now treated in
  -- a different way. For reference, I took a picture, 20201030_122839_DRO.jpg.
  --
  -- (-1) because the root height has an additive change.
  let jacobianL = Exp $ fromIntegral (nL - 1) * log xiL
      jacobianR = Exp $ fromIntegral (nR - 1) * log xiR
  return (tr', q, jacobianL * jacobianR)
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
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (HeightTree b)
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
  HeightTree a ->
  HeightTree a
scaleUltrametricTreeF h xi (Node _ lb ts) =
  Node () (lb & nodeHeightL .~ h') $ map (second $ nodeHeightL *~ xi') ts
  where
    xi' = either (error . (<>) "scaleSubTreeF:xi: ") id $ toHeight xi
    h' = either (error . (<>) "scaleSubTreeF:h: ") id $ toHeight h
