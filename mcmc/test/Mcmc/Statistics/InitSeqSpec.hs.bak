{- |
   Module      :  Mcmc.Statistics.InitSeqSpec
   Description :  Unit tests for Mcmc.Statistics.InitSeq
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

os2 :: Vector Double
os2 = V.fromList [2.39604875, 0.59133253, 0.49601200, 0.77828088, 0.36678584,
                  0.81875424, 0.59677332, 2.32285757, 0.13094589, 0.45258170,
                  0.10324180, 2.56827353, 0.50313726, 0.31220387, 1.33664320,
                  1.24568427, 0.80642147, 1.69867356, 0.03597119, 0.82510016,
                  0.95890956, 2.90166793, 0.67492229, 0.21564363]

xs2 :: Vector Int
xs2 = V.fromList [0,1,2,4,8,10,18,23]

ys2 :: Vector Double
ys2 = V.fromList [2.39604875,0.59133253,0.496012,0.36678584,0.13094589,0.1032418,3.597119e-2,0.21564363]

spec :: Spec
spec = do
  describe "gcm" $ do
    it "returns the greatest convex minorant" $ do
      gcm os `shouldBe` (xs, ys)
      gcm os2 `shouldBe` (xs2, ys2)
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
