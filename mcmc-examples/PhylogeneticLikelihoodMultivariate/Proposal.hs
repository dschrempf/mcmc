{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Proposal
-- Description :  Proposals on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Jul 23 09:10:07 2020.
module Proposal
  ( slideNode,
  )
where

import Control.Lens
import Mcmc.Proposal
import Mcmc.Proposal.Slide
import ELynx.Data.Tree
import System.Random.MWC

-- TODO: Think about how a (truncated) normal distribution could be used.

-- TODO: Think about how tuning could be enabled?

-- Minimum branch length.
eps :: Double
eps = 1e-8

-- Lens to a specific node.
nodeAt :: [Int] -> Lens' (Tree e a) (Tree e a)
nodeAt pth =
  lens
    (current . unsafeGoPath pth . fromTree)
    (\t t' -> let pos = unsafeGoPath pth . fromTree t in toTree $ pos {current = t'})

-- Lens to the branch of the root node.
rootBranch :: Lens' (Tree e a) e
rootBranch = lens branch (\(Node _ lb ts) br -> Node br lb ts)

modifyBranch :: (e -> e) -> Tree e a -> Tree e a
modifyBranch f (Node br lb ts) = Node (f br) lb ts

slideRootSample ::
  Tree Double a ->
  GenIO ->
  IO (Tree Double a)
slideRootSample (Node _ _ []) _ = error "slideRootSample: Cannot slide leaf node."
slideRootSample (Node br lb ts) g = do
  let br' = minimum $ map branch ts
  dx <- uniformR (negate $ br - eps, br' - eps) g
  return $ Node (br + dx) lb (map (modifyBranch (subtract dx)) ts)

-- slideNodeSample :: [Int] -> Tree Double a -> GenIO -> IO (Tree Double a)
-- slideNodeSample pth t g = case goPath pth $ fromTree t of
--   Nothing -> error $ "slideNodeSample: Could not find node with path " ++ show pth ++ "."
--   Just pos -> do
--     let ct = current pos
--     ct' <- slideRootSample ct g
--     return $ toTree $ insertTree ct' pos

slideRootSimple :: ProposalSimple (Tree Double a)
slideRootSimple = ProposalSimple slideRootSample Nothing

-- | Slide the node up and down using a uniform distribution truncated at
-- the parent node and the closest daughter node.
--
-- The node is specified by a path.
slideNode ::
  -- | Path to node on tree.
  [Int] ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  Proposal (Tree Double a)
slideNode pth n w = nodeAt pth >>> Proposal n w slideRootSimple Nothing

-- | Scale the branch of the node.
--
-- The node is specified by a path.
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
slideBranch pth n w s t = (nodeAt pth ^. rootBranch) >>> slideSymmetric n w s t
