-- |
-- Module      :  Normal
-- Description :  Benchmark Metropolis-Hastings-Green algorithm
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May  6 00:10:11 2020.
module Normal
  ( normalSlide,
    normalBactrianBench,
    normalLargeCycleBench,
    normalMC3,
  )
where

import Control.Monad
import Mcmc
import System.Random.MWC

trueMean :: Double
trueMean = 5

stdDev :: Double
stdDev = 4

lh :: LikelihoodFunction Double
lh = normal trueMean stdDev

cc :: Cycle Double
cc = cycleFromList [slideSymmetric 1.0 (PName "Medium") (PWeight 1) Tune]

mons :: [MonitorParameter Double]
mons = [monitorDouble "mu"]

monStd :: MonitorStdOut Double
monStd = monitorStdOut mons 200

mon :: Monitor Double
mon = Monitor monStd [] []

normalSlideBench :: GenIO -> IO ()
normalSlideBench g = do
  let s =
        Settings
          (AnalysisName "Normal")
          (BurnInWithAutoTuning 2000 500)
          (Iterations 20000)
          Overwrite
          Sequential
          NoSave
          Quiet
  a <- mhg noPrior lh cc mon 0 g
  void $ mcmc s a

ccLarge :: Cycle Double
ccLarge =
  cycleFromList
    [slideSymmetric 1.0 (PName $ "Medium " ++ show i) (PWeight 1) Tune | i <- [0 .. 100 :: Int]]

-- Should have the same run time as 'normalSlide'.
normalLargeCycleBench :: GenIO -> IO ()
normalLargeCycleBench g = do
  let s =
        Settings
          (AnalysisName "Normal")
          (BurnInWithAutoTuning 20 5)
          (Iterations 200)
          Overwrite
          Sequential
          NoSave
          Quiet
  a <- mhg noPrior lh ccLarge mon 0 g
  void $ mcmc s a

ccBactrian :: Cycle Double
ccBactrian = cycleFromList [slideBactrian 0.5 1.0 (PName "Bactrian") (PWeight 1) Tune]

normalBactrianBench :: GenIO -> IO ()
normalBactrianBench g = do
  let s =
        Settings
          (AnalysisName "NormalBactrian")
          (BurnInWithAutoTuning 2000 200)
          (Iterations 20000)
          Overwrite
          Sequential
          NoSave
          Quiet
  a <- mhg noPrior lh ccBactrian mon 0 g
  void $ mcmc s a

normalMC3 :: GenIO -> Int -> IO ()
normalMC3 g n = do
  let mcmcS =
        Settings
          (AnalysisName "MC3")
          (BurnInWithAutoTuning 200 20)
          (Iterations 2000)
          Overwrite
          Sequential
          NoSave
          Quiet
      mc3S = MC3Settings (NChains n) (SwapPeriod 2) (NSwaps 1)
  a <- mc3 mc3S noPrior lh cc mon 0 g
  void $ mcmc mcmcS a
