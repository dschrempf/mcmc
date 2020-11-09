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
-- ultrametric trees use a tree object directly storing the node height. For
-- this reason, __do not apply__ a mixed collection of __ultrametric__ and
-- __unconstrained__ proposals on the same tree object!
module Mcmc.Tree.Proposal.Ultrametric
  ( slideNodeUltrametric,
    scaleTreeUltrametric,
    scaleSubTreeUltrametric,
    pulleyUltrametric,
  )
where

import Control.Lens
import Control.Monad
import ELynx.Tree hiding (description)
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Tree.Lens
import Mcmc.Tree.Proposal.Common
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)
import Statistics.Distribution.Gamma
import System.Random.MWC

-- The branch is elongated by u. So if u is positive, the node height is
-- reduced.
slideNodeUltrametricF :: HasHeight a => Double -> Tree Length a -> Tree Length a
slideNodeUltrametricF u (Node br lb ts) =
  Node
    (br & lengthL %~ (+ u))
    (lb & heightL . lengthL %~ subtract u)
    (map (branchL . lengthL %~ subtract u) ts)

slideNodeUltrametricSimple ::
  HasHeight a =>
  -- Standard deviation.
  Double ->
  -- Tuning parameter.
  Double ->
  ProposalSimple (Tree Length a)
slideNodeUltrametricSimple _ _ (Node _ _ []) _ =
  error "slideNodeUltrametricSample: Cannot slide leaf node."
slideNodeUltrametricSimple s t tr@(Node br _ ts) g
  | br <= 0 =
    error $
      "slideNodeUltrametricSimple: Parent branch length is zero or negative: " ++ show br ++ "."
  | br' <= 0 =
    error $
      "slideNodeUltrametricSimple: Minimum branch length is zero or negative: " ++ show br' ++ "."
  | otherwise = do
    -- The determinant of the Jacobian is -1.0.
    (u, q) <- truncatedNormalSample s t a b g
    return (slideNodeUltrametricF u tr, q, 1.0)
  where
    br' = minimum $ map branch ts
    a = fromLength $ negate br
    b = fromLength br'

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
slideNodeUltrametric ::
  HasHeight a =>
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Tree Length a)
slideNodeUltrametric ds = createProposal description (slideNodeUltrametricSimple ds)
  where
    description = PDescription $ "Slide node ultrametric; sd: " ++ show ds

-- There is a distinction between tree types storing the node height only, and
-- ultrametric trees, because the Jacobian matrices differ. For ultrametric
-- trees, some branch lengths are constrained, whereas for normal trees, all
-- branch lengths can vary freely.
--
-- For the calculation of the Jacobian matrices of ultrametric trees, the height
-- of the inner nodes (and possibly the stem length) have to be used as
-- parameters.

-- TODO: Calculation of number of inner nodes may be slow.
nInnerNodes :: Tree e a -> Int
nInnerNodes (Node _ _ []) = 0
nInnerNodes tr = 1 + sum (map nInnerNodes $ forest tr)

scaleTreeUltrametricFunction :: HasHeight a => HandleStem -> Tree Length a -> Double -> Tree Length a
scaleTreeUltrametricFunction WithStem tr u =
  bimap
    (lengthL %~ (* u))
    (heightL . lengthL %~ (* u))
    tr
scaleTreeUltrametricFunction WithoutStem (Node br lb ts) u =
  Node
    br
    (lb & heightL . lengthL %~ (* u))
    $ map
      ( bimap
          (lengthL %~ (* u))
          (heightL . lengthL %~ (* u))
      )
      ts

-- (-2) for the entry corresponding to f(u)=1/u, (+1) for the stem makes
-- (-1) in total.
scaleTreeUltrametricJacobian :: Int -> HandleStem -> Tree e a -> Double -> Log Double
scaleTreeUltrametricJacobian n WithStem _ u = Exp $ fromIntegral (n - 1) * log u
scaleTreeUltrametricJacobian n WithoutStem _ u = Exp $ fromIntegral (n - 2) * log u

scaleTreeUltrametricSimple ::
  HasHeight a =>
  -- Number of inner nodes.
  Int ->
  HandleStem ->
  Double ->
  Double ->
  ProposalSimple (Tree Length a)
scaleTreeUltrametricSimple n s k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (scaleTreeUltrametricFunction s)
    (Just recip)
    (Just $ scaleTreeUltrametricJacobian n s)

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
--
-- __The height is changed.__ Do not use this proposal on a sub tree of an
-- ultrametric tree. Instead, use 'scaleSubTreeUltrametric'.
scaleTreeUltrametric ::
  HasHeight a =>
  -- | The tree is used to precompute the number of inner nodes for
  -- computational efficiency.
  Tree e b ->
  -- | Handle the stem?
  HandleStem ->
  -- | Shape.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Tree Length a)
scaleTreeUltrametric tr s k = createProposal description (scaleTreeUltrametricSimple n s k)
  where
    description = PDescription $ "Scale tree ultrametric; shape: " ++ show k
    n = nInnerNodes tr

-- The stem is elongated by u. So if u is positive, the node height is reduced.
--
-- The scaling factor for inner node heights is also given, since it is
-- calculated below.
slideBranchScaleSubTreeF :: HasHeight a => Double -> Double -> Tree Length a -> Tree Length a
slideBranchScaleSubTreeF u xi (Node br lb ts) =
  Node
    (br & lengthL %~ (+ u))
    (lb & heightL . lengthL %~ subtract u)
    $ map (bimap (lengthL %~ (* xi)) (heightL . lengthL %~ (* xi))) ts

scaleSubTreeUltrametricSimple ::
  (HasHeight a, Show a) =>
  -- Number of inner nodes.
  Int ->
  Double ->
  Double ->
  ProposalSimple (Tree Length a)
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
    ht = getHeight $ label tr
    a = fromLength $ negate br
    b = fromLength ht

-- | Scale the branches of the sub tree and slide the stem so that the tree
-- height is conserved.
--
-- See also 'scaleTreeUltrametric'.
--
-- A normal distribution truncated at the parent node (or the origin) and the
-- leaves is used to slide the given node.
scaleSubTreeUltrametric ::
  (HasHeight a, Show a) =>
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
scaleSubTreeUltrametric tr sd = createProposal description (scaleSubTreeUltrametricSimple n sd)
  where
    description = PDescription $ "Scale subtree ultrametrc; sd: " ++ show sd
    n = nInnerNodes tr

-- See 'pulleyTruncatedNormalSample'. However, we have to honor more constraints
-- in the ultrametric case.
pulleyUltrametricTruncatedNormalSample ::
  HasHeight a =>
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
    ht = getHeight lb
    -- The constraints are larger than 0.
    constraintRightBoundary = ht - brL
    constraintLeftBoundary = ht - brR
    a = fromLength $ negate $ minimum [brL, constraintLeftBoundary]
    b = fromLength $ minimum [brR, constraintRightBoundary]
pulleyUltrametricTruncatedNormalSample _ _ _ =
  error "pulleyUltrametricTruncatedNormalSample: Node is not bifurcating."

pulleyUltrametricSimple ::
  HasHeight a =>
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
  let hL = getHeight $ label l
      -- Scaling factor left. (hL - u)/hL = (1.0 - u/hL).
      xiL = 1.0 - u / fromLength hL
  -- Right.
  let hR = getHeight $ label r
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
  HasHeight a =>
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
