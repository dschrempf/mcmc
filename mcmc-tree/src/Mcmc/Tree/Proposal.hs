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
  ( -- * Slide branches
    slideBranch,
    slideNodeUltrametric,

    -- * Scale trees
    scaleTree,
    scaleTreeUltrametric,
    scaleSubTreeUltrametric,
    pulley,
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
import Mcmc.Proposal.Slide
import Mcmc.Tree.Height
import Mcmc.Tree.Lens
import Numeric.Log hiding (sum)
import Statistics.Distribution
import Statistics.Distribution.Gamma
import Statistics.Distribution.TruncatedNormal
import System.Random.MWC

-- | Slide branch.
--
-- Use a normal distribution with mean 0 and given standard deviation.
--
-- This proposal slides the root branch. To slide other branches, see 'subTreeAt'.
-- For example, @subTreeAt path @~ slideNodeUltrametric ...@.
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

-- The branch is elongated by u. So if u is positive, the node height is
-- reduced.
slideNodeUltrametricF :: HasHeight a => Double -> Tree Double a -> Tree Double a
slideNodeUltrametricF u (Node br lb ts) =
  Node (br + u) (applyHeight (subtract u) lb) (map (applyStem (subtract u)) ts)

slideNodeUltrametricSimple ::
  HasHeight a =>
  -- Standard deviation.
  Double ->
  -- Tuning parameter.
  Double ->
  ProposalSimple (Tree Double a)
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
    a = negate br
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
slideNodeUltrametric ::
  HasHeight a =>
  -- | Standard deviation.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
slideNodeUltrametric ds = createProposal description (slideNodeUltrametricSimple ds)
  where
    description = "Slide node ultrametric; sd: " ++ show ds

-- -- TODO.
-- data HandleStem = WithStem | WithoutStem
--
-- scaleTreeSimple :: HandleStem -> Double -> Double -> ProposalSimple (Tree Double a)

scaleTreeSimple :: Double -> Double -> ProposalSimple (Tree Double a)
scaleTreeSimple k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (\tr u -> first (* u) tr)
    (Just recip)
    (Just jac)
  where
    -- TODO: Scaling of the stem is included. This may be an issue.
    --
    -- TODO: Length calculation may be slow.
    jac tr u = Exp $ fromIntegral (length tr - 2) * log u

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
--
-- Because of the specificly used determinant of the Jacobian matrix, this
-- proposal is only valid, if all branch lengths free variables with strictly
-- positive values (including the stem). For example, ultrametric trees do not
-- fulfill this criterion.
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

-- TODO: There is now a distinction between tree types storing the node height
-- only, and ultrametric trees, because the Jacobian matrices differ. For
-- ultrametric trees, some branch lengths are constrained, whereas for normal
-- trees, all branch lengths can vary freely.
--
-- For the calculation of the Jacobian matrices of ultrametric trees, the height
-- of the inner nodes (and possibly the stem length) have to be used as
-- parameters.

-- TODO: Calculation of number of inner nodes may be slow.
nInnerNodes :: Tree e a -> Int
nInnerNodes (Node _ _ []) = 0
nInnerNodes tr = 1 + sum (map nInnerNodes $ forest tr)

scaleTreeUltrametricSimple ::
  HasHeight a =>
  Double ->
  Double ->
  ProposalSimple (Tree Double a)
scaleTreeUltrametricSimple k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (\tr x -> bimap (* x) (applyHeight (* x)) tr)
    (Just recip)
    (Just jac)
  where
    -- TODO: Scaling of the stem is included. This may be an issue.
    --
    -- (-2) for the entry corresponding to f(u)=1/u, (+1) for the stem makes
    -- (-1) in total.
    jac tr u = Exp $ fromIntegral (nInnerNodes tr - 1) * log u

-- | Scale all branches with a gamma distributed kernel of given shape. The
-- scale is set such that the mean is 1.0.
--
-- __The height is changed.__ Do not use this proposal on a sub tree of an
-- ultrametric tree. Instead, use 'scaleSubTreeUltrametric'.
scaleTreeUltrametric ::
  HasHeight a =>
  -- | Shape.
  Double ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
scaleTreeUltrametric k = createProposal description (scaleTreeUltrametricSimple k)
  where
    description = "Scale tree ultrametric; shape: " ++ show k

-- The stem is elongated by u. So if u is positive, the node height is reduced.
--
-- The scaling factor for inner node heights is also given, since it is
-- calculated below.
slideBranchScaleSubTreeF :: HasHeight a => Double -> Double -> Tree Double a -> Tree Double a
slideBranchScaleSubTreeF u xi (Node br lb ts) =
    Node (br + u) (applyHeight (subtract u) lb) $ map (bimap (* xi) (applyHeight (* xi))) ts

scaleSubTreeUltrametricSimple ::
  (HasHeight a, Show a) =>
  Double ->
  Double ->
  ProposalSimple (Tree Double a)
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
  let xi = 1.0 - u / ht
      jac = Exp $ fromIntegral (nInnerNodes tr - 1) * log xi
  return (slideBranchScaleSubTreeF u xi tr, q, jac)
  where
    br = branch tr
    ht = getHeight $ label tr
    a = negate br
    b = ht

-- | Scale the branches of the sub tree and slide the root branch so that the
-- tree height is conserved.
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
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
scaleSubTreeUltrametric sd = createProposal description (scaleSubTreeUltrametricSimple sd)
  where
    description = "Scale subtree ultrametrc; sd: " ++ show sd

-- XXX: I removed this function for now, because it has too many special cases.
-- The Jacobian differs for ultrametric and trees without constraints, and so we
-- would need at least three functions (UM, UM), (N, N), (N, U).

-- contra :: (Tree Double a, Tree Double b) -> Double -> (Tree Double a, Tree Double b)
-- contra (s, t) x = (first (* x) s, first (/ x) t)

-- scaleTreesContrarilySimple :: Double -> Double -> ProposalSimple (Tree Double a, Tree Double b)
-- scaleTreesContrarilySimple k t =
--   genericContinuous
--     (gammaDistr (k / t) (t / k))
--     contra
--     (Just recip)
--     (Just jac)
--   where
--     jac (l, r) u = undefined

-- -- | Scale all branches with a gamma distributed kernel of given shape. The
-- -- scale is set such that the mean is 1.0.
-- --
-- -- The two trees are scaled contrarily so that the product of their heights
-- -- stays constant. Contrary proposals are useful when parameters are confounded.
-- scaleTreesContrarily ::
--   -- | Shape.
--   Double ->
--   -- | Name.
--   String ->
--   -- | Weight.
--   Int ->
--   -- | Enable tuning.
--   Bool ->
--   Proposal (Tree Double a, Tree Double b)
-- scaleTreesContrarily k = createProposal description (scaleTreesContrarilySimple k)
--   where
--     description = "Scale trees contrarily; shape: " ++ show k

-- See 'truncatedNormalSample'. U is added to the left branch. I.e., if u is
-- positive, the left branch is elongated.
pulleyTruncatedNormalSample ::
  Double -> Double -> Tree Double a -> GenIO -> IO (Double, Log Double)
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
    a = negate brL
    b = brR
pulleyTruncatedNormalSample _ _ _ = error "pulleyTruncatedNormalSample: Node is not bifurcating."

pulleySimple :: Double -> Double -> ProposalSimple (Tree Double a)
pulleySimple s t tr@(Node br lb [l, r]) g = do
  (u, q) <- pulleyTruncatedNormalSample s t tr g
  let tr' = Node br lb [applyStem (+ u) l, applyStem (subtract u) r]
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
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
pulley s = createProposal description (pulleySimple s)
  where
    description = "Pulley; sd: " ++ show s

-- See 'pulleyTruncatedNormalSample'. However, we have to honor more constraints
-- in the ultrametric case.
pulleyUltrametricTruncatedNormalSample ::
  HasHeight a =>
  Double ->
  Double ->
  Tree Double a ->
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
    a = negate $ minimum [brL, constraintLeftBoundary]
    b = minimum [brR, constraintRightBoundary]
pulleyUltrametricTruncatedNormalSample _ _ _ =
  error "pulleyUltrametricTruncatedNormalSample: Node is not bifurcating."

pulleyUltrametricSimple :: HasHeight a => Double -> Double -> ProposalSimple (Tree Double a)
pulleyUltrametricSimple s t tr@(Node br lb [l, r]) g = do
  (u, q) <- pulleyUltrametricTruncatedNormalSample s t tr g
  -- Left.
  let hL = getHeight $ label l
      -- Scaling factor left. (hL - u)/hL = (1.0 - u/hL).
      xiL = 1.0 - u / hL
  -- Right.
  let hR = getHeight $ label r
      -- Scaling factor right. (hR + u)/hR = (1.0 + u/hR).
      xiR = 1.0 + u / hR
  let tr' = Node br lb [slideBranchScaleSubTreeF u xiL l, slideBranchScaleSubTreeF (negate u) xiR r]
  -- The calculation of the Jacobian matrix is very lengthy. Similar to before,
  -- we parameterize the tree into the two branch lengths leading to the pulley
  -- node, and the node heights of all other internal nodes. However, the left
  -- and right node heights are now treated in a different way. For reference, I
  -- took a picture, 20201030_122839_DRO.jpg.
  --
  -- (-1) because the root height has an additive change.
  let jacL = Exp $ fromIntegral (nInnerNodes l - 1) * log xiL
      jacR = Exp $ fromIntegral (nInnerNodes r - 1) * log xiR
  return (tr', q, jacL * jacR)
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
  String ->
  -- | Weight.
  Int ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
pulleyUltrametric d = createProposal description (pulleyUltrametricSimple d)
  where
    description = "Pulley ultrametric; sd: " ++ show d
