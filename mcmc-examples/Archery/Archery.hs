-- |
-- Module      :  Main
-- Description :  Modeling shots of an archer shots on a circular target
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 27 17:11:28 2020.
--
-- Please refer to the corresponding [RevBayes
-- tutorial](https://revbayes.github.io/tutorials/mcmc/archery.html).
module Main
  ( main,
  )
where

import Control.Monad
import Mcmc
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Exponential
import Statistics.Distribution.Gamma
import System.Random.MWC

-- State space, the accuracy of the archer.
type I = Double

-- Number of arrows.
nArrows :: Int
nArrows = 100

-- True precision of archer.
muTrue :: Double
muTrue = 1.0

-- Parameter of prior distribution.
alpha :: Double
alpha = 1.0

-- Likelihood function.
meanDistribution :: I -> GammaDistribution
meanDistribution x = gammaDistr (fromIntegral nArrows) (x / fromIntegral nArrows)

-- Simulated mean distance from center.
arrowMean :: GenIO -> IO I
arrowMean = genContVar (meanDistribution muTrue)

priorDistribution :: ExponentialDistribution
priorDistribution = exponential alpha

pr :: I -> Log Double
pr x
  | x <= 0 = pzero
  | otherwise = Exp $ logDensity priorDistribution x

-- Likelihood function.
lh :: I -> I -> Log Double
lh mu x
  | x <= 0 = pzero
  | otherwise = Exp $ logDensity (meanDistribution mu) x

moveCycle :: Cycle I
moveCycle = fromList [slideUniform "mu; slide uniform" 1 id 1.0 True]

monMu :: MonitorParameter I
monMu = monitorRealFloat "Mu" id

monStd :: MonitorStdOut I
monStd = monitorStdOut [monMu] 5000

monFile :: MonitorFile I
monFile = monitorFile "Mu" [monMu] 500

monMuBatch :: MonitorParameterBatch I
monMuBatch = monitorBatchMeanRealFloat "Mean mu" id

monBatch :: MonitorBatch I
monBatch = monitorBatch "Mu" [monMuBatch] 1000

mon :: Monitor I
mon = Monitor monStd [monFile] [monBatch]

nBurnIn :: Maybe Int
nBurnIn = Just 200000

nAutoTune :: Maybe Int
nAutoTune = Just 10000

nIter :: Int
nIter = 1000000

main :: IO ()
main = do
  g <- create
  mu <- arrowMean g
  putStrLn $ "True parameter: " <> show mu
  let s = noSave $ status "Archery" pr (lh mu) moveCycle mon 0.01 nBurnIn nAutoTune nIter g
  void $ mh s
