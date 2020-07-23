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

import Mcmc.Proposal
import ELynx.Data.Tree
import System.Random.MWC

-- TODO: Think about how a (truncated) normal distribution could be used.

-- TODO: Think about how tuning could be enabled?

-- Minimum branch length.
eps :: Double
eps = 1e-8

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

slideNodeSample :: [Int] -> Tree Double a -> GenIO -> IO (Tree Double a)
slideNodeSample pth t g = case goPath pth $ fromTree t of
  Nothing -> error $ "slideNodeSample: Could not find node with path " ++ show pth ++ "."
  Just pos -> do
    let ct = current pos
    ct' <- slideRootSample ct g
    return $ toTree $ insertTree ct' pos

slideNodeSimple :: [Int] -> ProposalSimple (Tree Double a)
slideNodeSimple pth = ProposalSimple (slideNodeSample pth) Nothing

-- | Slide the root node up and down using a uniform distribution truncated at
-- the origin and the closest daughter node.
slideNode ::
  -- | Path to node on tree.
  [Int] ->
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  Proposal (Tree Double a)
slideNode pth n w = Proposal n w (slideNodeSimple pth) Nothing
