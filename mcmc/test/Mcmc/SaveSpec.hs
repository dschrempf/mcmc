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
import Mcmc.Chain.Save
import Mcmc.Algorithm.Metropolis
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Directory
import qualified System.Random.MWC as R
import Test.Hspec

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

lh :: LikelihoodFunction Double
lh = Exp . logDensity (normalDistr trueMean trueStdDev)

proposals :: Cycle Double
proposals =
  cycleFromList
    [ slideSymmetric 0.1 (PName "Small") (PWeight 5) Tune,
      slideSymmetric 1.0 (PName "Medium") (PWeight 2) Tune,
      slideSymmetric 5.0 (PName "Large") (PWeight 2) Tune,
      slide 1.0 4.0 (PName "Skewed") (PWeight 1) Tune
    ]

monStd :: MonitorStdOut Double
monStd = monitorStdOut [monitorDouble "mu"] 10

mon :: Monitor Double
mon = Monitor monStd [] []

spec :: Spec
spec =
  describe "save and load" $
    it "doesn't change the MCMC chain" $
      do
        gen <- R.create
        let s = Settings "SaveSpec" (BurnInWithAutoTuning 20 10) 200 Overwrite NoSave Quiet
            c = fromMHG $ mhg noPrior lh proposals mon 0 gen
            c' = fromSavedChain noPrior lh proposals mon $ toSavedChain 100 c
        removeFile "SaveSpec.chain"
        r <- fromMHG <$> mcmc s (MHG c)
        r' <- fromMHG <$> mcmc s (MHG c')
        link r `shouldBe` link r'
        iteration r `shouldBe` iteration r'
        trace r `shouldBe` trace r'
        g <- R.save $ generator r
        g' <- R.save $ generator r'
        g `shouldBe` g'

-- -- TODO.
-- describe "mhContinue"
--   $ it "mh 200 + mhContinue 200 == mh 400"
--   $ do
--     gen1 <- create
--     let s1 = chain "SaveSpec" (const 1) likelihood proposals mon 0 nBurn nAutoTune 400 gen1
--     r1 <- mh s1
--     gen2 <- create
--     let s2 = chain "SaveSpec" (const 1) likelihood proposals mon 0 nBurn nAutoTune 200 gen2
--     r2' <- mh s2
--     r2 <- mhContinue 200 r2'
--     link r1 `shouldBe` link r2
--     iteration r1 `shouldBe` iteration r2
--     trace r1 `shouldBe` trace r2
--     g <- save $ generator r1
--     g' <- save $ generator r2
--     g `shouldBe` g'
