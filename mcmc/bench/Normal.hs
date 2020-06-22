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
nBurn = Just 2000

nAutoTune :: Maybe Int
nAutoTune = Just 200

nIter :: Int
nIter = 20000

normalBench :: GenIO -> IO ()
normalBench g = do
  let s = noSave $ status "Normal" (const 1) likelihood moveCycle mon 0 nBurn nAutoTune nIter g
  void $ mh s
