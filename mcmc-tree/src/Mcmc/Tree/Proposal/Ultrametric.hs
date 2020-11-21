-- |
-- Module      :  Mcmc.Tree.Proposal.Ultrametric
-- Description :  Proposals preserving tree height
-- Copyright   :  (c) Dominik Schrempf, 2020
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
module Mcmc.Tree.Proposal.Ultrametric
  ( -- TODO: slideAllNodesUltrametric
    slideNodeAtUltrametric,
    -- TODO: slideAllSubTreesUltrametric
    scaleSubTreeAtUltrametric,
    pulleyUltrametric,
  )
where

import Control.Lens
import Data.Bifunctor
import ELynx.Tree
import Mcmc.Proposal
import Mcmc.Tree.Lens
import Mcmc.Tree.Proposal.Common
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)
import System.Random.MWC

slideNodeAtUltrametricSimple ::
  Path ->
  -- Standard deviation.
  Double ->
  -- Tuning parameter.
  Double ->
  ProposalSimple HeightTree
slideNodeAtUltrametricSimple pth s t tr g
  | null ts = error $ "slideNodeAtUltrametricSimple: Cannot slide leaf: " ++ show nm ++ "."
  | otherwise = do
    -- The absolute value of the determinant of the Jacobian is 1.0.
    (h', q) <- truncatedNormalSample hNode s t hDaughter hParent g
    -- Use unsafe conversion.
    let setNodeHeight x = x & labelL . nodeHeightL .~ toLengthUnsafe h'
    return (toTree $ modifyTree setNodeHeight trPos, q, 1.0)
  where
    trPos = goPathUnsafe pth $ fromTree tr
    focus = current trPos
    ts = forest focus
    nm = nodeName $ label focus
    hNode = fromLength $ nodeHeight $ label focus
    hDaughter = fromLength $ maximum $ map (nodeHeight . label) ts
    hParent = fromLength $ nodeHeight $ label $ current $ goParentUnsafe trPos

-- | Slide node (for ultrametric trees).
--
-- For ultrametric trees, we cannot exclusively scale the branch such as with
-- 'Mcmc.Tree.Proposal.Unconstrained.scaleBranch', because this would change the
-- height and if the proposal is used on a non-root node, it would break
-- ultrametricity of the tree. Instead, we can slide the root node. That is,
-- when the stem is elongated, we need to shorten the daughter branches, and
-- vice versa, such that the tree height is conserved.
--
-- A normal distribution truncated at the origin and the closest daughter node
-- is used.
slideNodeAtUltrametric ::
  -- | Path to node.
  Path ->
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | Weight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal HeightTree
slideNodeAtUltrametric pth ds = createProposal description (slideNodeAtUltrametricSimple pth ds)
  where
    description = PDescription $ "Slide node ultrametric; sd: " ++ show ds

-- Calculate the number of inner nodes.
nInnerNodes :: Tree e a -> Int
nInnerNodes (Node _ _ []) = 0
nInnerNodes tr = 1 + sum (map nInnerNodes $ forest tr)

-- scaleTreeUltrametricFunction :: Measurable a => HandleStem -> Tree Length a -> Double -> Tree Length a
-- scaleTreeUltrametricFunction WithStem tr u =
--   bimap
--     (lengthL %~ (* u))
--     (measurableL . lengthL %~ (* u))
--     tr
-- scaleTreeUltrametricFunction WithoutStem (Node br lb ts) u =
--   Node
--     br
--     (lb & measurableL . lengthL %~ (* u))
--     $ map
--       ( bimap
--           (lengthL %~ (* u))
--           (measurableL . lengthL %~ (* u))
--       )
--       ts

-- -- (-2) for the entry corresponding to f(u)=1/u, (+1) for the stem makes
-- -- (-1) in total.
-- scaleTreeUltrametricJacobian :: Int -> HandleStem -> Tree e a -> Double -> Log Double
-- scaleTreeUltrametricJacobian n WithStem _ u = Exp $ fromIntegral (n - 1) * log u
-- scaleTreeUltrametricJacobian n WithoutStem _ u = Exp $ fromIntegral (n - 2) * log u

-- scaleTreeUltrametricSimple ::
--   Measurable a =>
--   -- Number of inner nodes.
--   Int ->
--   HandleStem ->
--   Double ->
--   Double ->
--   ProposalSimple (Tree Length a)
-- scaleTreeUltrametricSimple n s k t =
--   genericContinuous
--     (gammaDistr (k / t) (t / k))
--     (scaleTreeUltrametricFunction s)
--     (Just recip)
--     (Just $ scaleTreeUltrametricJacobian n s)

-- -- | Scale all branches with a gamma distributed kernel of given shape. The
-- -- scale is set such that the mean is 1.0.
-- --
-- -- __The height is changed.__ Do not use this proposal on a sub tree of an
-- -- ultrametric tree. Instead, use 'scaleSubTreeUltrametric'.
-- scaleTreeUltrametric ::
--   Measurable a =>
--   -- | The tree is used to precompute the number of inner nodes for
--   -- computational efficiency.
--   Tree e b ->
--   -- | Handle the stem?
--   HandleStem ->
--   -- | Shape.
--   Double ->
--   -- | Name.
--   PName ->
--   -- | Weight.
--   Weight ->
--   -- | Enable tuning.
--   Tune ->
--   Proposal (Tree Length a)
-- scaleTreeUltrametric tr s k = createProposal description (scaleTreeUltrametricSimple n s k)
--   where
--     description = PDescription $ "Scale tree ultrametric; shape: " ++ show k
--     n = nInnerNodes tr

-- The scaling factor for inner node heights is also given, since it is
-- calculated below.
scaleSubTreeF :: Double -> Double -> HeightTree -> HeightTree
scaleSubTreeF h xi (Node _ lb ts) =
  Node () (lb & nodeHeightL .~ h') $ map (second $ nodeHeightL *~ xi') ts
  where
    xi' = toLengthUnsafe xi
    h' = toLengthUnsafe h

scaleSubTreeAtUltrametricSimple ::
  -- Number of inner nodes.
  Int ->
  -- Path to sub tree.
  Path ->
  -- Standard deviation.
  Double ->
  -- Tuning parameter.
  Double ->
  ProposalSimple HeightTree
scaleSubTreeAtUltrametricSimple n pth ds t tr g
  | null ts =
    error $
      "scaleSubTreeAtUltrametricSimple: Cannot scale sub tree of leaf: " ++ show nm ++ "."
  | otherwise = do
    -- The determinant of the Jacobian is not included.
    (hNode', q) <- truncatedNormalSample hNode ds t 0 hParent g
    -- Scaling factor (xi, not x_i).
    let xi = hNode' / hNode
        -- (-1) because the root height has an additive change.
        jacobian = Exp $ fromIntegral (n - 1) * log xi
    return (scaleSubTreeF hNode' xi tr, q, jacobian)
  where
    trPos = goPathUnsafe pth $ fromTree tr
    focus = current trPos
    ts = forest focus
    nm = nodeName $ label focus
    hNode = fromLength $ nodeHeight $ label focus
    hParent = fromLength $ nodeHeight $ label $ current $ goParentUnsafe trPos

-- | Scale the node heights of the sub tree at given path.
--
-- A normal distribution truncated at the height of the parent node and the
-- leaves is used to determine the new height of the sub tree.
--
-- Call 'error' if
--
-- - this proposal is applied to the root node
--
-- - this proposal is applied to a leaf node
scaleSubTreeAtUltrametric ::
  -- | The tree is used to precompute the number of inner nodes for
  -- computational efficiency.
  Tree e b ->
  -- | Path to sub tree.
  Path ->
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | Weight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal HeightTree
scaleSubTreeAtUltrametric tr pth sd =
  createProposal
    description
    (scaleSubTreeAtUltrametricSimple n pth sd)
  where
    description = PDescription $ "Scale subtree ultrametrc; sd: " ++ show sd
    n = nInnerNodes $ current $ goPathUnsafe pth $ fromTree tr

-- TODO: This can be shortened.
--
-- See 'pulleyTruncatedNormalSample'. However, we have to honor more constraints
-- in the ultrametric case.
pulleyUltrametricTruncatedNormalSample ::
  Double ->
  Double ->
  HeightTree ->
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
    a = fromLength $ negate $ minimum [brL, constraintLeftBoundary]
    b = fromLength $ minimum [brR, constraintRightBoundary]
pulleyUltrametricTruncatedNormalSample _ _ _ =
  error "pulleyUltrametricTruncatedNormalSample: Node is not bifurcating."

pulleyUltrametricSimple ::
  -- Number of inner nodes of left tree.
  Int ->
  -- Number of inner nodes of right tree.
  Int ->
  Double ->
  Double ->
  ProposalSimple HeightTree
pulleyUltrametricSimple nL nR s t tr@(Node br lb [l, r]) g = do
  (u, q) <- pulleyUltrametricTruncatedNormalSample s t tr g
  -- Left.
  let hL = nodeHeight $ label l
      -- Scaling factor left. (hL - u)/hL = (1.0 - u/hL).
      xiL = 1.0 - u / fromLength hL
  -- Right.
  let hR = nodeHeight $ label r
      -- Scaling factor right. (hR + u)/hR = (1.0 + u/hR).
      xiR = 1.0 + u / fromLength hR
  let tr' = Node br lb [scaleSubTreeF u xiL l, scaleSubTreeF (negate u) xiR r]
  -- The calculation of the Jacobian matrix is very lengthy. Similar to before,
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
-- See 'Mcmc.Tree.Proposal.Unconstrained.pulley', but for ultrametric trees. The
-- sub trees are scaled so that the tree remains ultrametric.
pulleyUltrametric ::
  -- | The tree is used to precompute the number of inner nodes for
  -- computational efficiency.
  Tree e b ->
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | Weight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal HeightTree
pulleyUltrametric (Node _ _ [l, r]) d
  | null (forest l) = error "pulleyUltrametric: Left tree is a leaf."
  | null (forest r) = error "pulleyUltrametric: Right tree is a leaf."
  | otherwise = createProposal description (pulleyUltrametricSimple nL nR d)
  where
    description = PDescription $ "Pulley ultrametric; sd: " ++ show d
    nL = nInnerNodes l
    nR = nInnerNodes r
pulleyUltrametric _ _ = error "pulleyUltrametric: Node is not bifurcating."
