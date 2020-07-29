-- |
-- Module      :  Normal
-- Description :  Benchmark Metropolis-Hastings algorithm
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May  6 00:10:11 2020.
module Normal
  ( normalBench,
    normalBactrianBench,
  )
where

import Control.Monad
import Mcmc
import Numeric.Log as L
import Statistics.Distribution hiding
  ( mean,
    stdDev,
  )
import Statistics.Distribution.Normal
import System.Random.MWC

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

lh :: Double -> Log Double
lh = Exp . logDensity (normalDistr trueMean trueStdDev)

proposals :: Cycle Double
proposals =
  fromList
    [slideSymmetric "medium" 1 1.0 True]

mons :: [MonitorParameter Double]
mons = [monitorRealFloat "mu"]

monStd :: MonitorStdOut Double
monStd = monitorStdOut mons 200

-- monFile :: MonitorFile Double
-- monFile = monitorFile "Mu" mons 200

mon :: Monitor Double
mon = Monitor monStd [] []

nBurn :: Maybe Int
nBurn = Just 2000

nAutoTune :: Maybe Int
nAutoTune = Just 200

nIter :: Int
nIter = 20000

normalBench :: GenIO -> IO ()
normalBench g = do
  let s = quiet $ noSave $ status "Normal" (const 1) lh proposals mon 0 nBurn nAutoTune nIter g
  void $ mh s

proposalsBactrian :: Cycle Double
proposalsBactrian =
  fromList
    [slideBactrian "bactrian" 1 0.5 1.0 True]

normalBactrianBench :: GenIO -> IO ()
normalBactrianBench g = do
  let s = quiet $ noSave $ status "NormalBactrian" (const 1) lh proposalsBactrian mon 0 nBurn nAutoTune nIter g
  void $ mh s
