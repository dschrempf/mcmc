{- |
   Module      :  Mcmc.Statistics.InitSeqSpec
   Description :  Unit tests for Mcmc.Statistics.InitSeqSpec
   Copyright   :  (c) Dominik Schrempf, 2020
   License     :  GPL-3.0-or-later

   Maintainer  :  dominik.schrempf@gmail.com
   Stability   :  unstable
   Portability :  portable

Creation date: Fri Jun  5 17:01:22 2020.

-}

module Mcmc.Statistics.InitSeqSpec
  ( spec
  ) where

import           Test.Hspec
import qualified Data.Vector.Unboxed as V
import Data.Vector.Unboxed (Vector)

import Mcmc.Statistics.InitSeq

os :: Vector Double
os = V.fromList [19, 16, 17, 12, 19, 6, 13, 1, 4, 0]

xs :: Vector Int
xs = V.fromList [0, 1, 7, 9]

ys :: Vector Double
ys = V.fromList [19.0, 16.0, 1.0, 0.0]

ysSmooth :: Vector Double
ysSmooth = V.fromList [19.0, 16.0, 13.5, 11.0, 8.5, 6.0, 3.5, 1.0, 0.5, 0.0]

xsWeird :: Vector Int
xsWeird = V.fromList [-2001, 7293, 9128]

ysWeird :: Vector Double
ysWeird = V.fromList [88, 22, 920938]

spec :: Spec
spec = do
  describe "gcm" $ do
    it "returns the greatest convex minorant" $
      gcm os `shouldBe` (xs, ys)
    it "shouldn't fail on a singleton vector" $
      gcm (V.singleton 1.0) `shouldBe` (V.singleton 0, V.singleton 1.0)
    it "shouldn't fail on an empty vector" $
      gcm V.empty `shouldBe` (V.empty, V.empty)
  describe "smooth" $ do
    it "correctly smooths various test cases" $ do
      smooth xs ys `shouldBe` ysSmooth
      let val = smooth xsWeird ysWeird
      V.length val `shouldBe` 9128 + 2001 + 1
      val V.! 0 `shouldBe` 88
      val V.! (7293 + 2001) `shouldBe` 22
      V.last val `shouldBe` 920938
    it "shouldn't fail on a singleton vector" $
      smooth (V.singleton 0) (V.singleton 1.0) `shouldBe` V.singleton 1.0
    it "shouldn't fail on an empty vector" $
      smooth V.empty V.empty `shouldBe` V.empty
