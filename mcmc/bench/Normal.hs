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
    [slideSymmetric 1.0 (PName "Medium") (PWeight 1) Tune]

mons :: [MonitorParameter Double]
mons = [monitorDouble "mu"]

monStd :: MonitorStdOut Double
monStd = monitorStdOut mons 200

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
  let
    e = quiet def
    c = chain "Normal" (const 1) lh proposals mon 0 nBurn nAutoTune nIter g
  void $ mh e c

proposalsBactrian :: Cycle Double
proposalsBactrian =
  fromList
    [slideBactrian 0.5 1.0 (PName "Bactrian") (PWeight 1) Tune]

normalBactrianBench :: GenIO -> IO ()
normalBactrianBench g = do
  let e = quiet def
      c = chain "NormalBactrian" (const 1) lh proposalsBactrian mon 0 nBurn nAutoTune nIter g
  void $ mh e c
