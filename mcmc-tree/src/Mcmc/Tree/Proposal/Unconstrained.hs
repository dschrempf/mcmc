-- |
-- Module      :  Mcmc.Tree.Proposal.Unconstrained
-- Description :  Proposals on trees with unconstrained branch lengths
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 09:10:07 2020.
module Mcmc.Tree.Proposal.Unconstrained
  ( scaleBranch,
    scaleTree,
    pulley,
  )
where

import Control.Lens
import Data.Bifunctor
import ELynx.Tree hiding (description)
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Proposal.Scale
import Mcmc.Tree.Lens
import Mcmc.Tree.Types
import Mcmc.Tree.Proposal.Common
import Numeric.Log hiding (sum)
import Statistics.Distribution.Gamma
import System.Random.MWC

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

