-- |
-- Module      :  Tools
-- Description :  Tools for MCMC runs on trees
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri May 28 11:57:07 2021.
module Tools
  ( scaleBothBranches,
  )
where

import ELynx.Tree
import Mcmc
import Mcmc.Proposal
import Statistics.Distribution
import Statistics.Distribution.Gamma

mult' :: Double -> Length -> Length
mult' u = toLengthUnsafe . (* u) . fromLength

scaleBothBranchesSimple :: Shape -> TuningParameter -> ProposalSimple (Tree Length a)
scaleBothBranchesSimple k t (Node cb cl [Node lb ll ls, Node rb rl rs]) g = do
  u <- genContVar d g
  let x' = Node cb cl [Node (mult' u lb) ll ls, Node (mult' u rb) rl rs]
      u1 = recip u
      qXY = Exp $ logDensity d u
      qYX = Exp $ logDensity d u1
      j = Exp $ log u1
  return (x', qYX / qXY, j)
  where
    d = gammaDistr (k / t) (recip k * t)
scaleBothBranchesSimple _ _ _ _ = error "scaleBothBranchesSimple: Tree is not bifurcating."

-- | Scale both daughter branches of the root.
--
-- Assume those branches are constrained (not two free variables).
scaleBothBranches :: Shape -> PName -> PWeight -> Tune -> Proposal (Tree Length a)
scaleBothBranches k = createProposal dscr (scaleBothBranchesSimple k) dim
  where
    dscr = PDescription "Scale daugther branches of root"
    dim = PDimension 1
