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
import Mcmc.Chain.Chain
import Numeric.Log
import Statistics.Distribution hiding
  ( mean,
    stdDev,
  )
import Statistics.Distribution.Normal
import qualified System.Random.MWC as R
import Test.Hspec

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

lh :: Double -> Log Double
lh = Exp . logDensity (normalDistr trueMean trueStdDev)

proposals :: Cycle Double
proposals =
  fromList
    [ slideSymmetric 0.1 (PName "Small") (PWeight 5) Tune,
      slideSymmetric 1.0 (PName "Medium") (PWeight 2) Tune,
      slideSymmetric 5.0 (PName "Large") (PWeight 2) Tune,
      slide 1.0 4.0 (PName "Skewed") (PWeight 1) Tune
    ]

monStd :: MonitorStdOut Double
monStd = monitorStdOut [monitorDouble "mu"] 10

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
  describe "save and load" $
    it "doesn't change the MCMC chain" $
      do
        gen <- R.create
        let env = Environment Force (SaveN 100) Quiet
            ch = chain "SaveSpec" (const 1) lh proposals mon 0 nBurn nAutoTune nIter gen
        save env ch
        (env', ch') <- load (const 1) lh proposals mon Nothing "SaveSpec"
        r <- mh env ch
        r' <- mh env' ch'
        -- Done during 'loadStatus'.
        -- removeFile "SaveSpec.json"
        item r `shouldBe` item r'
        iteration r `shouldBe` iteration r'
        trace r `shouldBe` trace r'
        g <- R.save $ generator r
        g' <- R.save $ generator r'
        g `shouldBe` g'
        env `shouldBe` env'

-- -- TODO: Splitmix. This will only work with a splittable generator
-- -- because getNIterations changes the generator.
-- describe "mhContinue"
--   $ it "mh 200 + mhContinue 200 == mh 400"
--   $ do
--     gen1 <- create
--     let s1 = status "SaveSpec" (const 1) likelihood proposals mon 0 nBurn nAutoTune 400 gen1
--     r1 <- mh s1
--     gen2 <- create
--     let s2 = status "SaveSpec" (const 1) likelihood proposals mon 0 nBurn nAutoTune 200 gen2
--     r2' <- mh s2
--     r2 <- mhContinue 200 r2'
--     item r1 `shouldBe` item r2
--     iteration r1 `shouldBe` iteration r2
--     trace r1 `shouldBe` trace r2
--     g <- save $ generator r1
--     g' <- save $ generator r2
--     g `shouldBe` g'
