-- |
-- Module      :  Normal
-- Description :  Benchmark Metropolis-Hastings-Green algorithm
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May  6 00:10:11 2020.
module Normal
  ( normalSlideBench,
    normalBactrianBench,
    normalLargeCycleBench,
    normalMC3,
  )
where

import Control.Monad
import Mcmc
import System.Random.Stateful

trueMean :: Double
trueMean = 5

stdDev :: Double
stdDev = 4

lh :: LikelihoodFunction Double
lh = normal trueMean stdDev

cc :: Cycle Double
cc = cycleFromList [slideSymmetric 1.0 (PName "Medium") (pWeight 1) Tune]

mons :: [MonitorParameter Double]
mons = [monitorDouble "mu"]

monStd :: MonitorStdOut Double
monStd = monitorStdOut mons 200

mon :: Monitor Double
mon = Monitor monStd [] []

normalSlideBench :: IOGenM StdGen -> IO ()
normalSlideBench g = do
  let s =
        Settings
          (AnalysisName "Normal")
          (BurnInWithAutoTuning 2000 500)
          (Iterations 20000)
          TraceAuto
          Overwrite
          Sequential
          NoSave
          LogStdOutOnly
          Quiet
  a <- mhg s noPrior lh cc mon 0 g
  void $ mcmc s a

ccLarge :: Cycle Double
ccLarge =
  cycleFromList
    [slideSymmetric 1.0 (PName $ "Medium " ++ show i) (pWeight 1) Tune | i <- [0 .. 100 :: Int]]

-- Should have the same run time as 'normalSlide'.
normalLargeCycleBench :: IOGenM StdGen -> IO ()
normalLargeCycleBench g = do
  let s =
        Settings
          (AnalysisName "Normal")
          (BurnInWithAutoTuning 20 5)
          (Iterations 200)
          TraceAuto
          Overwrite
          Sequential
          NoSave
          LogStdOutOnly
          Quiet
  a <- mhg s noPrior lh ccLarge mon 0 g
  void $ mcmc s a

ccBactrian :: Cycle Double
ccBactrian = cycleFromList [slideBactrian 0.5 1.0 (PName "Bactrian") (pWeight 1) Tune]

normalBactrianBench :: IOGenM StdGen -> IO ()
normalBactrianBench g = do
  let s =
        Settings
          (AnalysisName "NormalBactrian")
          (BurnInWithAutoTuning 2000 200)
          (Iterations 20000)
          TraceAuto
          Overwrite
          Sequential
          NoSave
          LogStdOutOnly
          Quiet
  a <- mhg s noPrior lh ccBactrian mon 0 g
  void $ mcmc s a

normalMC3 :: IOGenM StdGen -> Int -> IO ()
normalMC3 g n = do
  let mcmcS =
        Settings
          (AnalysisName "MC3")
          (BurnInWithAutoTuning 200 20)
          (Iterations 2000)
          TraceAuto
          Overwrite
          Sequential
          NoSave
          LogStdOutOnly
          Quiet
      mc3S = MC3Settings (NChains n) (SwapPeriod 2) (NSwaps 1)
  a <- mc3 mc3S mcmcS noPrior lh cc mon 0 g
  void $ mcmc mcmcS a
