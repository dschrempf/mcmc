-- |
--   Module      :  Mcmc.ProposalSpec
--   Description :  Unit tests for Mcmc.Proposal
--   Copyright   :  2021 Dominik Schrempf
--   License     :  GPL-3.0-or-later
--
--   Maintainer  :  dominik.schrempf@gmail.com
--   Stability   :  unstable
--   Portability :  portable
--
-- Creation date: Thu Jun 25 11:46:05 2020.
module Mcmc.ProposalSpec
  ( spec,
  )
where

import Mcmc.Cycle
import Mcmc.Proposal
import Mcmc.Proposal.Slide
import System.Random.MWC
import Test.Hspec

p1 :: Proposal Double
p1 = slideSymmetric 1.0 (PName "Test 1") (pWeight 1) Tune

p2 :: Proposal Double
p2 = slideSymmetric 1.0 (PName "Test 2") (pWeight 3) Tune

c :: Cycle Double
c = cycleFromList [p1, p2]

spec :: Spec
spec =
  describe "prepareProposals" $
    it "returns the correct number of proposals in a cycle" $
      do
        g <- create
        l1 <- length <$> prepareProposals AllProposals c g
        l1 `shouldBe` 4
        l2 <- length <$> prepareProposals AllProposals (setOrder RandomReversibleO c) g
        l2 `shouldBe` 8
        o3 <- prepareProposals AllProposals (setOrder SequentialReversibleO c) g
        length o3 `shouldBe` 8
        o3 == [p1, p2, p2, p2, p2, p2, p2, p1] `shouldBe` True
