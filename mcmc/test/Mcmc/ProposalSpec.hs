-- |
--   Module      :  Mcmc.ProposalSpec
--   Description :  Unit tests for Mcmc.Proposal
--   Copyright   :  (c) Dominik Schrempf, 2020
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

import Mcmc.Proposal
import Mcmc.Proposal.Slide
import System.Random.MWC
import Test.Hspec

p1 :: Proposal Double
p1 = slideSymmetric "test1" 1 1.0 True

p2 :: Proposal Double
p2 = slideSymmetric "test2" 3 1.0 True

c :: Cycle Double
c = fromList [p1, p2]

spec :: Spec
spec =
  describe "getNIterations" $
    it "returns the correct number of proposals in a cycle" $
      do
        g <- create
        l1 <- length . head <$> getNIterations c 1 g
        l1 `shouldBe` 4
        l2 <- length . head <$> getNIterations (setOrder RandomReversibleO c) 1 g
        l2 `shouldBe` 8
        o3 <- head <$> getNIterations (setOrder SequentialReversibleO c) 1 g
        length o3 `shouldBe` 8
        o3 `shouldBe` [p1, p2, p2, p2, p2, p2, p2, p1]
