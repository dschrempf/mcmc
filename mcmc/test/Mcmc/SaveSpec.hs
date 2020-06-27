-- |
--   Module      :  Mcmc.SaveSpec
--   Description :  Unit tests for Mcmc.Save
--   Copyright   :  (c) Dominik Schrempf, 2020
--   License     :  GPL-3.0-or-later
--
--   Maintainer  :  dominik.schrempf@gmail.com
--   Stability   :  unstable
--   Portability :  portable
--
-- Creation date: Tue Jun 16 14:32:55 2020.
module Mcmc.SaveSpec
  ( spec,
  )
where

import Mcmc
import Mcmc.Save
import Mcmc.Status hiding (save)
import Numeric.Log
import Statistics.Distribution hiding
  ( mean,
    stdDev,
  )
import Statistics.Distribution.Normal
import System.Directory
import System.Random.MWC
import Test.Hspec

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

lh :: Double -> Log Double
lh = Exp . logDensity (normalDistr trueMean trueStdDev)

moveCycle :: Cycle Double
moveCycle =
  fromList
    [ slideSymmetric "small" 5 id 0.1 True,
      slideSymmetric "medium" 2 id 1.0 True,
      slideSymmetric "large" 2 id 5.0 True,
      slide "skewed" 1 id 1.0 4.0 True
    ]

monStd :: MonitorStdOut Double
monStd = monitorStdOut [monitorRealFloat "mu" id] 10

mon :: Monitor Double
mon = Monitor monStd [] []

nBurn :: Maybe Int
nBurn = Just 20

nAutoTune :: Maybe Int
nAutoTune = Just 10

nIter :: Int
nIter = 200

spec :: Spec
spec =
  describe "saveStatus and loadStatus"
    $ it "doesn't change the MCMC chain"
    $ do
      gen <- create
      let s = quiet $ noSave $ status "SaveSpec" (const 1) lh moveCycle mon 0 nBurn nAutoTune nIter gen
      saveStatus "SaveSpec.json" s
      s' <- loadStatus (const 1) lh moveCycle mon "SaveSpec.json"
      removeFile "SaveSpec.json"
      r <- mh s
      r' <- mh s'
      item r `shouldBe` item r'
      iteration r `shouldBe` iteration r'
      trace r `shouldBe` trace r'
      g <- save $ generator r
      g' <- save $ generator r'
      g `shouldBe` g'

-- -- TODO: Splitmix. This will only work with a splittable generator
-- -- because getNCycles changes the generator.
-- describe "mhContinue"
--   $ it "mh 200 + mhContinue 200 == mh 400"
--   $ do
--     gen1 <- create
--     let s1 = noSave $ status "SaveSpec" (const 1) likelihood moveCycle mon 0 nBurn nAutoTune 400 gen1
--     r1 <- mh s1
--     gen2 <- create
--     let s2 = noSave $ status "SaveSpec" (const 1) likelihood moveCycle mon 0 nBurn nAutoTune 200 gen2
--     r2' <- mh s2
--     r2 <- mhContinue 200 r2'
--     item r1 `shouldBe` item r2
--     iteration r1 `shouldBe` iteration r2
--     trace r1 `shouldBe` trace r2
--     g <- save $ generator r1
--     g' <- save $ generator r2
--     g `shouldBe` g'
