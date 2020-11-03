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
-- For reasons of computational efficiency some functions, for example, the ones
-- working with ultrametric trees, use a tree object directly storing the node height.
module Mcmc.Tree.Proposal
  ( -- * Unconstrained trees
    scaleBranch,
    scaleTree,
    pulley,

    -- * Ultrametric trees
    slideNodeUltrametric,
    scaleTreeUltrametric,
    scaleSubTreeUltrametric,
    pulleyUltrametric,
  )
where

import Control.Lens
import Control.Monad
import Data.Bifunctor
-- import Debug.Trace
import ELynx.Tree hiding (description)
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Proposal.Scale
import Mcmc.Tree.Lens
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)
import Statistics.Distribution
import Statistics.Distribution.Gamma
import Statistics.Distribution.TruncatedNormal
import System.Random.MWC

-- Which lens to use?
lengthL :: Lens' Length Double
-- -- For debugging, around ten percent slower.
-- lengthL = lengthE
-- For production, faster but does not fail on negative branch lengths.
lengthL = lengthU

-- | Scale branch.
--
-- See 'scaleUnbiased'.
--
-- This proposal scales the stem. To slide other branches, see 'subTreeAtE'. For
-- example, @subTreeAtE path @~ slideNodeUltrametric ...@.
scaleBranch ::
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Tree Length a)
scaleBranch s n w t = (branchL . lengthL) @~ scaleUnbiased s n w t

-- -- Minimum branch length.
-- eps :: Double
-- eps = 1e-12

-- A very specific function that samples a delta value from the truncated normal
-- distribution with given bounds [a,b] and also computes the required factor of
-- the Metropolis-Hastings proposal ratio.
--
-- NO JACOBIAN IS COMPUTED, because we do not know how the proposal will be used.
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
  u <- genContinuous d g
  -- Compute Metropolis-Hastings factor.
  let a' = a - u
      b' = b - u
      d' = either error id $ truncatedNormalDistr 0 s' a' b'
      qXY = Exp $ logDensity d u
      qYX = Exp $ logDensity d' (- u)
  -- NO JACOBIAN IS COMPUTED.
  return (u, qYX / qXY)

-- TODO: CONTINUE HERE.

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
-- 'scaleBranch', because this would change the height and if the proposal is
-- used on a non-root node, it would break ultrametricity of the tree. Instead,
-- we can slide the root node. That is, when the stem is elongated, we need to
-- shorten the daughter branches, and vice versa, such that the tree height is
-- conserved.
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

scaleTreeFunction :: HandleStem -> Tree Length a -> Double -> Tree Length a
scaleTreeFunction WithStem tr u = first (lengthL %~ (* u)) tr
scaleTreeFunction WithoutStem (Node br lb ts) u =
  Node br lb $ map (first (lengthL %~ (* u))) ts

-- TODO: Length calculation may be slow.
scaleTreeJacobian :: HandleStem -> Tree e a -> Double -> Log Double
scaleTreeJacobian WithStem tr u = Exp $ fromIntegral (length tr - 2) * log u
scaleTreeJacobian WithoutStem tr u = Exp $ fromIntegral (length tr - 3) * log u

scaleTreeSimple :: HandleStem -> Double -> Double -> ProposalSimple (Tree Length a)
scaleTreeSimple s k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (scaleTreeFunction s)
    (Just recip)
    (Just $ scaleTreeJacobian s)

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
--
-- Because of the specificly used determinant of the Jacobian matrix, this
-- proposal is only valid, if all branch lengths free variables with strictly
-- positive values (including the stem). For example, ultrametric trees do not
-- fulfill this criterion.
scaleTree ::
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
scaleTree s k = createProposal description (scaleTreeSimple s k)
  where
    description = PDescription $ "Scale tree; shape: " ++ show k

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
scaleTreeUltrametricJacobian :: HandleStem -> Tree e a -> Double -> Log Double
scaleTreeUltrametricJacobian WithStem tr u = Exp $ fromIntegral (nInnerNodes tr - 1) * log u
scaleTreeUltrametricJacobian WithoutStem tr u = Exp $ fromIntegral (nInnerNodes tr - 2) * log u

scaleTreeUltrametricSimple ::
  HasHeight a =>
  HandleStem ->
  Double ->
  Double ->
  ProposalSimple (Tree Length a)
scaleTreeUltrametricSimple s k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (scaleTreeUltrametricFunction s)
    (Just recip)
    (Just $ scaleTreeUltrametricJacobian s)

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
--
-- __The height is changed.__ Do not use this proposal on a sub tree of an
-- ultrametric tree. Instead, use 'scaleSubTreeUltrametric'.
scaleTreeUltrametric ::
  HasHeight a =>
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
scaleTreeUltrametric s k = createProposal description (scaleTreeUltrametricSimple s k)
  where
    description = PDescription $ "Scale tree ultrametric; shape: " ++ show k

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
  Double ->
  Double ->
  ProposalSimple (Tree Length a)
scaleSubTreeUltrametricSimple _ _ (Node _ _ []) _ =
  error "scaleSubTreeUltrametricSample: Cannot scale sub tree of leaf node."
scaleSubTreeUltrametricSimple ds t tr g = do
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
  -- Scaling factor (xi, not x_i) (ht - u)/ht = (1.0 - u/ht).
  let xi = 1.0 - u / b
      jacobian = Exp $ fromIntegral (nInnerNodes tr - 1) * log xi
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
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Tree Length a)
scaleSubTreeUltrametric sd = createProposal description (scaleSubTreeUltrametricSimple sd)
  where
    description = PDescription $ "Scale subtree ultrametrc; sd: " ++ show sd

-- See 'truncatedNormalSample'. U is added to the left branch. I.e., if u is
-- positive, the left branch is elongated.
pulleyTruncatedNormalSample ::
  Double -> Double -> Tree Length a -> GenIO -> IO (Double, Log Double)
pulleyTruncatedNormalSample s t (Node _ _ [l, r])
  | brL <= 0 =
    error $
      "pulleyTruncatedNormalSample: Left branch is zero or negative: " ++ show brL ++ "."
  | brR <= 0 =
    error $
      "pulleyTruncatedNormalSample: Right branch is zero or negative: " ++ show brR ++ "."
  | otherwise = truncatedNormalSample s t a b
  where
    brL = branch l
    brR = branch r
    a = fromLength $ negate brL
    b = fromLength brR
pulleyTruncatedNormalSample _ _ _ = error "pulleyTruncatedNormalSample: Node is not bifurcating."

pulleySimple :: Double -> Double -> ProposalSimple (Tree Length a)
pulleySimple s t tr@(Node br lb [l, r]) g = do
  (u, q) <- pulleyTruncatedNormalSample s t tr g
  let tr' =
        Node
          br
          lb
          [ l & branchL . lengthL %~ (+ u),
            r & branchL . lengthL %~ subtract u
          ]
  -- The determinant of the Jacobian matrix is (-1).
  return (tr', q, 1.0)
pulleySimple _ _ _ _ = error "pulleySimple: Node is not bifurcating."

-- | Use a node as a pulley.
--
-- For bifurcating nodes; change the daughter branch lengths contrarily. The
-- first and second daughter branches are elongated/shortened by the same
-- amount.
--
-- The heights of the two sub trees change.
pulley ::
  -- | Standard deviation.
  Double ->
  -- | Name.
  PName ->
  -- | PWeight.
  PWeight ->
  -- | Enable tuning.
  Tune ->
  Proposal (Tree Length a)
pulley s = createProposal description (pulleySimple s)
  where
    description = PDescription $ "Pulley; sd: " ++ show s

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
  | otherwise = do truncatedNormalSample s t a b
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

pulleyUltrametricSimple :: HasHeight a => Double -> Double -> ProposalSimple (Tree Length a)
pulleyUltrametricSimple s t tr@(Node br lb [l, r]) g = do
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
  let jacobianL = Exp $ fromIntegral (nInnerNodes l - 1) * log xiL
      jacobianR = Exp $ fromIntegral (nInnerNodes r - 1) * log xiR
  return (tr', q, jacobianL * jacobianR)
pulleyUltrametricSimple _ _ _ _ = error "pulleyUltrametricSimple: Node is not bifurcating."

-- | Use a node as a pulley.
--
-- See 'pulley', but for ultrametric trees. The sub trees are scaled such that
-- the tree heights are conserved and the tree remains ultrametric.
pulleyUltrametric ::
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
pulleyUltrametric d = createProposal description (pulleyUltrametricSimple d)
  where
    description = PDescription $ "Pulley ultrametric; sd: " ++ show d
