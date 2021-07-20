-- |
--   Module      :  Mcmc.SaveSpec
--   Description :  Unit tests for Mcmc.Save
--   Copyright   :  (c) Dominik Schrempf, 2021
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
import Mcmc.Chain.Trace
import Statistics.Distribution
import Statistics.Distribution.Normal
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
    [ slideSymmetric 0.1 (PName "Small") (pWeight 5) Tune,
      slideSymmetric 1.0 (PName "Medium") (pWeight 2) Tune,
      slideSymmetric 5.0 (PName "Large") (pWeight 2) Tune,
      slide 1.0 4.0 (PName "Skewed") (pWeight 1) Tune
    ]

monStd :: MonitorStdOut Double
monStd = monitorStdOut [monitorDouble "mu"] 10

mon :: Monitor Double
mon = Monitor monStd [] []

settings :: Settings
settings =
  Settings
    (AnalysisName "SaveSpec")
    (BurnInWithAutoTuning 20 10)
    (Iterations 100)
    TraceAuto
    Overwrite
    Sequential
    NoSave
    LogStdOutOnly
    Quiet

spec :: Spec
spec = do
  describe "save and load" $
    it "doesn't change the MCMC chain" $
      do
        gen <- R.create
        a <- mhg settings noPrior lh proposals mon 0 gen
        c <- fromMHG <$> mcmc settings a
        savedChain <- toSavedChain c
        c' <- fromSavedChain noPrior lh proposals mon savedChain
        putStrLn "@load . save@ should be @id@."
        link c `shouldBe` link c'
        iteration c `shouldBe` iteration c'
        frozenT1 <- freezeT (trace c)
        frozenT1' <- freezeT (trace c')
        frozenT1 `shouldBe` frozenT1'
        -- g1 <- R.save $ generator c
        -- g1' <- R.save $ generator c'
        -- g1 `shouldBe` g1'
        putStrLn "Sampling from the chains should be the same."
        r <- fromMHG <$> mcmcContinue (Iterations 100) settings (MHG c)
        r' <- fromMHG <$> mcmcContinue (Iterations 100) settings (MHG c')
        link r `shouldBe` link r'
        iteration r `shouldBe` iteration r'
        frozenT2 <- freezeT (trace c)
        frozenT2' <- freezeT (trace c')
        frozenT2 `shouldBe` frozenT2'
        g2 <- R.save $ generator r
        g2' <- R.save $ generator r'
        g2 `shouldBe` g2'

  describe "mhContinue" $
    it "mcmc 50 + mcmcContinue 50 == mcmc 100" $
      do
        gen1 <- R.create
        a1 <- mhg settings noPrior lh proposals mon 0 gen1
        r1 <- fromMHG <$> mcmc settings a1
        gen2 <- R.create
        let settings' = settings {sIterations = Iterations 50}
        a2 <- mhg settings' noPrior lh proposals mon 0 gen2
        r2' <- mcmc settings' a2
        r2 <- fromMHG <$> mcmcContinue (Iterations 50) settings' r2'
        link r1 `shouldBe` link r2
        iteration r1 `shouldBe` iteration r2
        frozenT1 <- freezeT (trace r1)
        frozenT2 <- freezeT (trace r2)
        frozenT1 `shouldBe` frozenT2
        g <- R.save $ generator r1
        g' <- R.save $ generator r2
        g `shouldBe` g'
