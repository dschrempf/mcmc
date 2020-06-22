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

import Mcmc hiding (save)
import Numeric.Log
import Statistics.Distribution hiding
  ( mean,
    stdDev,
  )
import Statistics.Distribution.Normal
import System.Random.MWC
import Test.Hspec

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

likelihood :: Double -> Log Double
likelihood = Exp . logDensity (normalDistr trueMean trueStdDev)

moveCycle :: Cycle Double
moveCycle =
  fromList
    [ slideSymmetric "small" 5 id 0.1 True,
      slideSymmetric "medium" 2 id 1.0 True,
      slideSymmetric "large" 2 id 5.0 True,
      slide "skewed" 1 id 1.0 4.0 True
    ]

monStd :: MonitorStdOut Double
monStd = monitorStdOut [monitorRealFloat "mu" id] 200

mon :: Monitor Double
mon = Monitor monStd [] []

nBurn :: Maybe Int
nBurn = Just 20

nAutoTune :: Maybe Int
nAutoTune = Just 10

nIter :: Int
nIter = 200

spec :: Spec
spec = describe "saveStatus and loadStatus"
  $ it "doesn't change the MCMC chain"
  $ do
    gen <- create
    let s = status "SaveSpec" (const 1) likelihood moveCycle mon 0 nBurn nAutoTune nIter gen
    saveStatus "SaveSpec.json" s
    s' <- loadStatus (const 1) likelihood moveCycle mon "SaveSpec.json"
    r <- mh s
    r' <- mh s'
    g <- save $ generator r
    g' <- save $ generator r'
    item r `shouldBe` item r'
    iteration r `shouldBe` iteration r'
    trace r `shouldBe` trace r'
    g `shouldBe` g'
