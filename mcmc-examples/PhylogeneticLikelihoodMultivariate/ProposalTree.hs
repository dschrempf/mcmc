{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  ProposalTree
-- Description :  Proposals on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 09:10:07 2020.
module ProposalTree
  ( slideNode,
    slideBranch,
    scaleTree,
  )
where

import Control.Monad.Primitive
import Data.Bifunctor
import Control.Lens
import ELynx.Data.Tree
import Mcmc.Proposal
import Mcmc.Proposal.Generic
import Mcmc.Proposal.Slide
import Numeric.Log
import Statistics.Distribution.Gamma
import System.Random.MWC
import System.Random.MWC.Distributions

-- Minimum branch length.
eps :: Double
eps = 1e-8

-- Lens to a specific node.
nodeAt :: [Int] -> Lens' (Tree e a) (Tree e a)
nodeAt pth =
  lens
    (current . unsafeGoPath pth . fromTree)
    ( \t t' ->
        let pos = unsafeGoPath pth $ fromTree t
         in toTree $ pos {current = t'}
    )

-- Lens to the branch of the root node.
rootBranch :: Lens' (Tree e a) e
rootBranch = lens branch (\(Node _ lb ts) br -> Node br lb ts)

-- Be careful, this will loop forever if the parameters are not well chosen.
truncatedNormal :: PrimMonad m => Double -> Double -> Double -> Double -> Gen (PrimState m) -> m Double
truncatedNormal a b m s g = do
  x' <- normal m s g
  case x' of
    x
      | x < a -> truncatedNormal a b m s g
      | x > b -> truncatedNormal a b m s g
      | otherwise -> return x

modifyBranch :: (e -> e) -> Tree e a -> Tree e a
modifyBranch f (Node br lb ts) = Node (f br) lb ts

slideRootSample ::
  Double ->
  Tree Double a ->
  GenIO ->
  IO (Tree Double a, Log Double)
slideRootSample _ (Node _ _ []) _ = error "slideRootSample: Cannot slide leaf node."
slideRootSample t (Node br lb ts) g = do
  let br' = minimum $ map branch ts
      a = negate $ br - eps
      b = br' - eps
      -- Don't let the standard deviation be too high, because then the normal
      -- variable will be rejected many times in 'truncatedNormal'.
      s = min (b - a) (t / 2 * (b - a))
  dx <- truncatedNormal a b 0 s g
  let t' = Node (br + dx) lb (map (modifyBranch (subtract dx)) ts)
  return (t', 1.0)

slideRootSimple :: Double -> ProposalSimple (Tree Double a)
slideRootSimple t = ProposalSimple $ slideRootSample t

-- | Slide the node up and down using a uniform distribution truncated at
-- the parent node and the closest daughter node.
--
-- The node to slide is specified by a path.
slideNode ::
  -- | Path to node on tree.
  [Int] ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Tune the move.
  Bool ->
  Proposal (Tree Double a)
slideNode pth n w t = nodeAt pth @~ createProposal n w slideRootSimple t

-- | Slide the branch of the node.
--
-- The branch to slide is specified by a path to the node.
slideBranch ::
  -- | Path to node on tree.
  [Int] ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Standard deviation.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
slideBranch pth n w s t = (nodeAt pth . rootBranch) @~ slideSymmetric n w s t

scaleTreeSimple :: Double -> Double -> Double -> ProposalSimple (Tree Double a)
scaleTreeSimple k th t =
  proposalGenericContinuous
    (gammaDistr (k / t) (th * t))
    (\tr x -> first (* x) tr)
    (Just recip)

-- | Scale all branches of a tree with a Gamma distributed kernel.
scaleTree ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Shape.
  Double ->
  -- | Scale.
  Double ->
  -- | Enable tuning.
  Bool ->
  Proposal (Tree Double a)
scaleTree n w k th = createProposal n w (scaleTreeSimple k th)
