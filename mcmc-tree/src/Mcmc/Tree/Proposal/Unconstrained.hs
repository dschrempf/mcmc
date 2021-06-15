-- |
-- Module      :  Mcmc.Tree.Proposal.Unconstrained
-- Description :  Proposals on trees with unconstrained branch lengths
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 09:10:07 2020.
--
-- For proposals on ultrametric trees, see "Mcmc.Tree.Proposal.Ultrametric".
module Mcmc.Tree.Proposal.Unconstrained
  ( scaleBranch,
    scaleBranches,
    scaleTree,
    scaleSubTrees,
    pulley,
    scaleNormAndTreeContrarily,
    scaleVarianceAndTree,

    -- * Helper functions
    scaleUnconstrainedTreeF,
    scaleUnconstrainedTreeWithoutStemF,
    scaleUnconstrainedStem,
  )
where

import Control.Lens
import Data.Bifunctor
import ELynx.Tree
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Proposal.Scale
import Mcmc.Statistics.Types
import Mcmc.Tree.Import ()
import Mcmc.Tree.Lens
import Mcmc.Tree.Proposal.Common
import Mcmc.Tree.Types
import Numeric.Log hiding (sum)
import Statistics.Distribution.Gamma
import System.Random.MWC

-- | Scale branch.
--
-- This proposal scales the stem of the given tree. To slide inner branches, use
-- 'subTreeAtL' or 'scaleBranches'. For example,
--
-- @subTreeAtL path @~ scaleBranch ...@.
--
-- See 'scaleUnbiased'.
scaleBranch ::
  Shape ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Tree Length a)
scaleBranch s n w t = (branchL . lengthL) @~ scaleUnbiased s n w t

-- | Scale the branches of a given tree.
--
-- See 'scaleBranch'.
scaleBranches ::
  Tree e a ->
  HandleLayer ->
  Shape ->
  -- | Base name of proposals.
  PName ->
  PWeight ->
  Tune ->
  [Proposal (Tree Length b)]
scaleBranches tr hd s n w t =
  [ subTreeAtL pth
      @~ scaleBranch s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      -- Filter layers.
      hd $ length pth
  ]
  where
    name lb = n <> PName (" branch " ++ show lb)

scaleTreeJacobian ::
  -- Number of branches.
  Int ->
  a ->
  Double ->
  Jacobian
scaleTreeJacobian n _ u = Exp $ fromIntegral (n - 2) * log u

scaleTreeSimple ::
  -- Number of branches.
  Int ->
  Shape ->
  TuningParameter ->
  ProposalSimple (Tree Length a)
scaleTreeSimple n k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (flip scaleUnconstrainedTreeF)
    (Just recip)
    (Just $ scaleTreeJacobian n)

-- | Scale all branches including the stem.
--
-- A gamma distributed kernel of given shape is used. The scale is set such that
-- the mean is 1.0.
--
-- NOTE: Because the determinant of the Jacobian matrix depends on the number of
-- branches scaled, this proposal is only valid if all branch lengths including
-- the stem are unconstrained and strictly positive.
--
-- Call 'error' if:
--
-- - A branch length is zero or negative.
scaleTree ::
  -- | The topology of the tree is used to precompute the number of inner nodes.
  Tree e b ->
  Shape ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Tree Length a)
scaleTree tr k = createProposal description (scaleTreeSimple n k) (PDimension n)
  where
    description = PDescription $ "Scale tree; shape: " ++ show k
    n = length tr

-- | Scale the sub trees of a given tree.
--
-- See 'scaleTree'.
--
-- The weights are assigned as described in
-- 'Mcmc.Tree.Proposal.Ultrametric.scaleSubTreesUltrametric'.
--
-- Do not scale the leaves.
scaleSubTrees ::
  Tree e a ->
  HandleLayer ->
  Shape ->
  -- | Base name of proposals.
  PName ->
  -- | Minimum weight.
  PWeight ->
  -- | Maximum weight.
  PWeight ->
  Tune ->
  [Proposal (Tree Length b)]
scaleSubTrees tr hd s n wMin wMax t =
  [ subTreeAtL pth
      @~ scaleTree focus s (name lb) w t
    | (pth, lb) <- itoList $ identify tr,
      let focus = tr ^. subTreeAtL pth,
      let currentLayer = length pth,
      let currentDepth = depth focus,
      -- Subtract 2 because leaves have depth one and are not scaled.
      let w = pWeight $ minimum [fromPWeight wMin + currentDepth - 2, fromPWeight wMax],
      -- Do not scale the leaves, because 'scaleBranch' is faster.
      not $ null $ forest focus,
      -- Filter other layers.
      hd currentLayer
  ]
  where
    name lb = n <> PName (" node " ++ show lb)

-- See 'truncatedNormalSample'. U is added to the left branch. I.e., if u is
-- positive, the left branch is elongated.
pulleyTruncatedNormalSample ::
  StandardDeviation -> TuningParameter -> Tree Length a -> GenIO -> IO (Double, Log Double)
pulleyTruncatedNormalSample s t (Node _ _ [l, r])
  | brL <= 0 =
    error $
      "pulleyTruncatedNormalSample: Left branch is zero or negative: " ++ show brL ++ "."
  | brR <= 0 =
    error $
      "pulleyTruncatedNormalSample: Right branch is zero or negative: " ++ show brR ++ "."
  | otherwise = truncatedNormalSample 0 s t a b
  where
    brL = branch l
    brR = branch r
    a = fromLength $ negate brL
    b = fromLength brR
pulleyTruncatedNormalSample _ _ _ = error "pulleyTruncatedNormalSample: Node is not bifurcating."

pulleySimple :: StandardDeviation -> TuningParameter -> ProposalSimple (Tree Length a)
pulleySimple s t tr@(Node br lb [l, r]) g = do
  (u, q) <- pulleyTruncatedNormalSample s t tr g
  let tr' =
        Node
          br
          lb
          [ l & branchL . lengthL +~ u,
            r & branchL . lengthL -~ u
          ]
  -- The determinant of the Jacobian matrix is (-1).
  return (tr', q, 1.0)
pulleySimple _ _ _ _ = error "pulleySimple: Node is not bifurcating."

-- | Use a node as a pulley.
--
-- Change the daughter branch lengths contrarily. The left and right daughter
-- branches are elongated/shortened by the same amount.
--
-- The heights of the two sub trees change.
--
-- Call 'error' if:
--
-- - The node is not bifurcating.
pulley ::
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Tree Length a)
pulley s = createProposal description (pulleySimple s) (PDimension 2)
  where
    description = PDescription $ "Pulley; sd: " ++ show s

scaleNormAndTreeContrarilyFunction ::
  (Double, Tree Length a) ->
  Double ->
  (Double, Tree Length a)
-- Do not touch the stem, because this leads to problems with negative stem
-- lengths (or NaNs).
scaleNormAndTreeContrarilyFunction (x, tr) u = (x / u, scaleUnconstrainedTreeWithoutStemF u tr)

scaleNormAndTreeContrarilySimple ::
  -- Number of branches.
  Int ->
  Shape ->
  TuningParameter ->
  ProposalSimple (Double, Tree Length a)
scaleNormAndTreeContrarilySimple n k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    scaleNormAndTreeContrarilyFunction
    (Just recip)
    (Just jacobianFunction)
  where
    -- Minus 2 because of reverse transform.
    -- Minus 1 because of scaling the norm contrarily.
    jacobianFunction _ u = Exp $ fromIntegral (n - 2 - 1) * log u

-- | Let the branch lengths of an unconstrained tree be distributed according to
-- a distribution with a given mean hyper-parameter. This proposal scales the
-- mean hyper-parameter and the branches.
--
-- NOTE: Because the determinant of the Jacobian matrix depends on the number of
-- branches scaled, this proposal is only valid if all inner branch lengths are
-- unconstrained and strictly positive. The stem is ignored.
--
-- See also 'scaleVarianceAndTree'.
scaleNormAndTreeContrarily ::
  -- | The topology of the tree is used to precompute the number of inner branches.
  Tree e a ->
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Double, Tree Length b)
scaleNormAndTreeContrarily tr sd =
  createProposal
    description
    (scaleNormAndTreeContrarilySimple nBranches sd)
    (PDimension $ nBranches + 1)
  where
    description = PDescription $ "Scale norm and tree contrarily; sd: " <> show sd
    -- Ignore the stem, see note above.
    nBranches = length tr - 1

scaleVarianceAndTreeFunction ::
  Int ->
  (Double, Tree Length a) ->
  Double ->
  (Double, Tree Length a)
-- Do not touch the stem, this leads to problems with negative stem lengths (or
-- NaNs).
scaleVarianceAndTreeFunction n (x, tr) u =
  (x * u * u, tr & forestL %~ map (first (lengthL %~ f)))
  where
    -- The calculation of sample mean requires a separate traversal of the tree.
    -- This may be slow.
    --
    -- Specifically ignore the stem.
    s = fromLength $ sum $ concatMap branches $ forest tr
    -- Sample mean. We assume that 'n' does not count the stem.
    mu = s / fromIntegral n
    -- Force NaN when new value is negative. This avoids numerical problems with
    -- negative branch lengths.
    f b = let b' = (b - mu) * u + mu in if b' > 0 then b' else 0 / 0

scaleVarianceAndTreeSimple ::
  -- Number of branches.
  Int ->
  Shape ->
  TuningParameter ->
  ProposalSimple (Double, Tree Length a)
scaleVarianceAndTreeSimple n k t =
  genericContinuous
    (gammaDistr (k / t) (t / k))
    (scaleVarianceAndTreeFunction n)
    (Just recip)
    (Just jacobianFunction)
  where
    n1 = recip (fromIntegral n)
    -- Minus 2 because of the reverse transform.
    -- Plus 2 because of the scaling of the variance.
    -- Is 0 :).
    -- However, n times a complicated factor because the mean is used.
    jacobianFunction _ u = Exp $ fromIntegral n * log (u - n1 * u + n1)

-- | Let the branch lengths of an unconstrained tree be distributed according to
-- a distribution with a given variance hyper-parameter. This proposal scales
-- the variance hyper-parameter and the branches. The branches are scaled such
-- that the sample variance is scaled in the same way as the variance
-- hyper-parameter.
--
-- NOTE: Because the determinant of the Jacobian matrix depends on the number of
-- branches scaled, this proposal is only valid if all inner branch lengths are
-- unconstrained and strictly positive. The stem is ignored.
--
-- NOTE: The sample mean is used to scale the parameters. This has two
-- disadvantages: (1) It is slow, and (2) the sample mean may be off
-- considerably. However, the actual mean of the distribution is, in general,
-- unknown.
scaleVarianceAndTree ::
  -- | The topology of the tree is used to precompute the number of inner branches.
  Tree e a ->
  StandardDeviation ->
  PName ->
  PWeight ->
  Tune ->
  Proposal (Double, Tree Length b)
scaleVarianceAndTree tr sd =
  createProposal
    description
    (scaleVarianceAndTreeSimple nBranches sd)
    (PDimension $ nBranches + 1)
  where
    description = PDescription $ "Scale variance and tree; sd: " <> show sd
    -- Ignore the stem, see note above.
    nBranches = length tr - 1

-- | Scale the branches of an unconstrained tree.
scaleUnconstrainedTreeF :: Double -> Tree Length a -> Tree Length a
scaleUnconstrainedTreeF u = first (lengthL *~ u)

-- | Scale the branches of an unconstrained tree. Do not scale the stem.
scaleUnconstrainedTreeWithoutStemF :: Double -> Tree Length a -> Tree Length a
scaleUnconstrainedTreeWithoutStemF u tr = tr & forestL %~ map (first $ lengthL *~ u)

-- | Scale the stem of an unconstrained tree.
scaleUnconstrainedStem :: Double -> Tree Length a -> Tree Length a
scaleUnconstrainedStem u tr = tr & branchL %~ lengthL *~ u
