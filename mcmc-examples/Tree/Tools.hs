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

scaleBothBranchesSimple :: TuningParameter -> ProposalSimple (Tree Length a)
scaleBothBranchesSimple t x g = undefined

-- | Scale both daughter branches of the root.
scaleBothBranches :: PName -> PWeight -> Tune -> Proposal (Tree Length a)
scaleBothBranches = createProposal dscr scaleBothBranchesSimple dim
  where dscr = PDescription "Scale daugther branches of root"
        dim = PDimension 1
