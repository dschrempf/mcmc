-- |
--   Module      :  Mcmc.MoveSpec
--   Description :  Unit tests for Mcmc.Move
--   Copyright   :  (c) Dominik Schrempf, 2020
--   License     :  GPL-3.0-or-later
--
--   Maintainer  :  dominik.schrempf@gmail.com
--   Stability   :  unstable
--   Portability :  portable
--
-- Creation date: Thu Jun 25 11:46:05 2020.
module Mcmc.MoveSpec
  ( spec,
  )
where

import Mcmc.Move
import Mcmc.Move.Slide
import System.Random.MWC
import Test.Hspec

mv1 :: Move Double
mv1 = slideSymmetric "test1" 1 id 1.0 True

mv2 :: Move Double
mv2 = slideSymmetric "test2" 3 id 1.0 True

c :: Cycle Double
c = fromList [mv1, mv2]

spec :: Spec
spec =
  describe "getNCycles"
    $ it "returns the correct number of moves in a cycle"
    $ do
      g <- create
      l1 <- length . head <$> getNCycles c 1 g
      l1 `shouldBe` 4
      l2 <- length . head <$> getNCycles (setOrder RandomReversibleO c) 1 g
      l2 `shouldBe` 8
      o3 <- head <$> getNCycles (setOrder SequentialReversibleO c) 1 g
      length o3 `shouldBe` 8
      o3 `shouldBe` [mv1, mv2, mv2, mv2, mv2, mv2, mv2, mv1]
