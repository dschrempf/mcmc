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

spec :: Spec
spec = do
  describe "save and load" $
    it "doesn't change the MCMC chain" $
      do
        gen <- R.create
        let s =
              Settings
                (AnalysisName "SaveSpec")
                (BurnInWithAutoTuning 20 10)
                (Iterations 200)
                Overwrite
                Sequential
                NoSave
                LogStdOutOnly
                Quiet
        c <- fromMHG <$> mhg noPrior lh proposals mon TraceAuto 0 gen
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
        r <- fromMHG <$> mcmc s (MHG c)
        r' <- fromMHG <$> mcmc s (MHG c')
        link r `shouldBe` link r'
        iteration r `shouldBe` iteration r'
        frozenT2 <- freezeT (trace c)
        frozenT2' <- freezeT (trace c')
        frozenT2 `shouldBe` frozenT2'
        g2 <- R.save $ generator r
        g2' <- R.save $ generator r'
        g2 `shouldBe` g2'

-- -- TODO: 'mhContinue'.
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
