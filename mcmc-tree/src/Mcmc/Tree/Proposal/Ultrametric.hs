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
  ( slideNodeAtUltrametric,
    -- scaleTreeUltrametric,
    scaleSubTreeUltrametric,
    pulleyUltrametric,
  )
where

import Control.Lens
import Control.Monad
import ELynx.Tree
import Mcmc.Proposal
import Mcmc.Tree.Lens
import Mcmc.Tree.Proposal.Common
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)
import System.Random.MWC

-- setNodeHeight :: Length -> HeightTree -> HeightTree
-- setNodeHeight h (Node () lb ts) = Node () (lb & nodeHeight .~ h) ts

slideNodeAtUltrametricSimple ::
  Path ->
  -- Standard deviation.
  Double ->
  -- Tuning parameter.
  Double ->
  ProposalSimple HeightTree
slideNodeAtUltrametricSimple _ _ _ (Node _ _ []) _ =
  error "slideNodeAtUltrametricSample: Cannot slide leaf node."
slideNodeAtUltrametricSimple [] _ _ (Node _ _ []) _ =
  error "slideNodeAtUltrametricSample: Cannot slide root node."
slideNodeAtUltrametricSimple pth s t tr g
  | a >= b =
    error $
      "slideNodeAtUltrametricSimple: Maximum height of daughters "
        ++ show a
        ++ " is larger equal parent height "
        ++ show b
        ++ "."
  | otherwise = do
    -- The determinant of the Jacobian is -1.0.
    (h', q) <- truncatedNormalSample s t a b g
    -- Use unsafe conversion.
    let setNodeHeight x = x & labelL . nodeHeight .~ toLengthUnsafe h'
    return (toTree $ modifyTree setNodeHeight trPos, q, 1.0)
  where
    trPos = goPathUnsafe pth $ fromTree tr
    -- The maximum of the node heights of the daughter nodes.
    a = fromLength $ maximum $ map (_nodeHeight . label) $ forest $ current trPos
    -- The parent node height.
    b = fromLength $ _nodeHeight $ label $ current $ goParentUnsafe trPos

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
  -- | PWeight.
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
--   -- | PWeight.
--   PWeight ->
--   -- | Enable tuning.
--   Tune ->
--   Proposal (Tree Length a)
-- scaleTreeUltrametric tr s k = createProposal description (scaleTreeUltrametricSimple n s k)
--   where
--     description = PDescription $ "Scale tree ultrametric; shape: " ++ show k
--     n = nInnerNodes tr

-- The stem is elongated by u. So if u is positive, the node height is reduced.
--
-- The scaling factor for inner node heights is also given, since it is
-- calculated below.
slideBranchScaleSubTreeF :: Double -> Double -> HeightTree -> HeightTree
slideBranchScaleSubTreeF u xi (Node br lb ts) =
  Node
    (br & lengthL %~ (+ u))
    (lb & measurableL . lengthL %~ subtract u)
    $ map (bimap (lengthL %~ (* xi)) (measurableL . lengthL %~ (* xi))) ts

scaleSubTreeUltrametricSimple ::
  -- Number of inner nodes.
  Int ->
  Double ->
  Double ->
  ProposalSimple HeightTree
scaleSubTreeUltrametricSimple _ _ _ (Node _ _ []) _ =
  error "scaleSubTreeUltrametricSample: Cannot scale sub tree of leaf node."
scaleSubTreeUltrametricSimple n ds t tr g = do
  when
    (br <= 0)
    ( error $
        "scaleSubTreeUltrametricSimple: Parent branch length is zero or negative: " ++ show br ++ "." ++ show tr
    )
  when
    (ht <= 0)
    ( error $
        "scaleSubTreeUltrametricSimple: Node height is zero or negative: " ++ show ht ++ "."
    )
  -- The determinant of the Jacobian is not included.
  (u, q) <- truncatedNormalSample ds t a b g
  -- For the calculation of the Jacobian matrix, parameterize the tree into the
  -- stem length, the height of the root, and the heights of other inner nodes.
  --
  -- (-1) because the root height has an additive change.
  --
  -- Scaling factor (xi, not x_i) is (ht - u)/ht = (1.0 - u/ht).
  let xi = 1.0 - u / b
      jacobian = Exp $ fromIntegral (n - 1) * log xi
  return (slideBranchScaleSubTreeF u xi tr, q, jacobian)
  where
    br = branch tr
    ht = getLen $ label tr
    a = fromLength $ negate br
    b = fromLength ht

-- | Scale the branches of the sub tree and slide the stem so that the tree
-- height is conserved.
--
-- A normal distribution truncated at the parent node (or the origin) and the
-- leaves is used to slide the given node.
scaleSubTreeUltrametric ::
  -- | The tree is used to precompute the number of inner nodes for
  -- computational efficiency.
  Tree e b ->
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal HeightTree
scaleSubTreeUltrametric tr sd = createProposal description (scaleSubTreeUltrametricSimple n sd)
  where
    description = PDescription $ "Scale subtree ultrametrc; sd: " ++ show sd
    n = nInnerNodes tr

-- See 'pulleyTruncatedNormalSample'. However, we have to honor more constraints
-- in the ultrametric case.
pulleyUltrametricTruncatedNormalSample ::
  Measurable a =>
  Double ->
  Double ->
  Tree Length a ->
  GenIO ->
  IO (Double, Log Double)
pulleyUltrametricTruncatedNormalSample s t (Node _ lb [l, r])
  | brL <= 0 =
    error $
      "pulleyUltrametricTruncatedNormalSample: Left branch is zero or negative: " ++ show brL ++ "."
  | brR <= 0 =
    error $
      "pulleyUltrametricTruncatedNormalSample: Right branch is zero or negative: " ++ show brR ++ "."
  | otherwise = truncatedNormalSample s t a b
  where
    -- Left and right branch length.
    brL = branch l
    brR = branch r
    -- The new branch lengths are not allowed to exceed the height of the node.
    ht = getLen lb
    -- The constraints are larger than 0.
    constraintRightBoundary = ht - brL
    constraintLeftBoundary = ht - brR
    a = fromLength $ negate $ minimum [brL, constraintLeftBoundary]
    b = fromLength $ minimum [brR, constraintRightBoundary]
pulleyUltrametricTruncatedNormalSample _ _ _ =
  error "pulleyUltrametricTruncatedNormalSample: Node is not bifurcating."

pulleyUltrametricSimple ::
  Measurable a =>
  -- Number of inner nodes of left tree.
  Int ->
  -- Number of inner nodes of right tree.
  Int ->
  Double ->
  Double ->
  ProposalSimple (Tree Length a)
pulleyUltrametricSimple nL nR s t tr@(Node br lb [l, r]) g = do
  (u, q) <- pulleyUltrametricTruncatedNormalSample s t tr g
  -- Left.
  let hL = getLen $ label l
      -- Scaling factor left. (hL - u)/hL = (1.0 - u/hL).
      xiL = 1.0 - u / fromLength hL
  -- Right.
  let hR = getLen $ label r
      -- Scaling factor right. (hR + u)/hR = (1.0 + u/hR).
      xiR = 1.0 + u / fromLength hR
  let tr' = Node br lb [slideBranchScaleSubTreeF u xiL l, slideBranchScaleSubTreeF (negate u) xiR r]
  -- The calculation of the Jacobian matrix is very lengthy. Similar to before,
  -- we parameterize the tree into the two branch lengths leading to the pulley
  -- node, and the node heights of all other internal nodes. However, the left
  -- and right node heights are now treated in a different way. For reference, I
  -- took a picture, 20201030_122839_DRO.jpg.
  --
  -- (-1) because the root height has an additive change.
  let jacobianL = Exp $ fromIntegral (nL - 1) * log xiL
      jacobianR = Exp $ fromIntegral (nR - 1) * log xiR
  return (tr', q, jacobianL * jacobianR)
pulleyUltrametricSimple _ _ _ _ _ _ = error "pulleyUltrametricSimple: Node is not bifurcating."

-- | Use a node as a pulley.
--
-- See 'Mcmc.Tree.Proposal.Unconstrained.pulley', but for ultrametric trees. The
-- sub trees are scaled such that the tree heights are conserved and the tree
-- remains ultrametric.
pulleyUltrametric ::
  Measurable a =>
  -- | The tree is used to precompute the number of inner nodes for
  -- computational efficiency.
  Tree e b ->
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Tree Length a)
pulleyUltrametric (Node _ _ [l, r]) d
  | null (forest l) = error "pulleyUltrametric: Left tree is a leaf."
  | null (forest r) = error "pulleyUltrametric: Right tree is a leaf."
  | otherwise = createProposal description (pulleyUltrametricSimple nL nR d)
  where
    description = PDescription $ "Pulley ultrametric; sd: " ++ show d
    nL = nInnerNodes l
    nR = nInnerNodes r
pulleyUltrametric _ _ = error "pulleyUltrametric: Node is not bifurcating."
