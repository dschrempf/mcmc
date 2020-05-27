{- |
Module      :  Metropolis
Description :  Benchmark Metropolis-Hastings algorithm
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 00:10:11 2020.

-}

module Metropolis
  ( mhBench
  ) where

import qualified Data.Vector.Unboxed as V
import Numeric.Log as L
import Statistics.Distribution hiding (mean, stdDev)
import Statistics.Distribution.Normal
import Statistics.Sample
import System.Random.MWC

import Statistics.Mcmc

import Statistics.Mcmc.Acceptance
import Statistics.Mcmc.Item
import Statistics.Mcmc.Trace

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

posterior :: Double -> Log Double
posterior = Exp . logDensity (normalDistr trueMean trueStdDev)

moveCycle :: Cycle Double
moveCycle = fromList [ (slideDouble "small"  0   0.1 True, 5)
                     , (slideDouble "medium" 0   1.0 True, 2)
                     , (slideDouble "large"  0   5.0 True, 2)
                     , (slideDouble "skewed" 1.0 4.0 True, 1)]

summarize :: [Double] -> (Double, Double)
summarize xs = (mean v, stdDev v)
  where v = V.fromList xs

nBurn :: Maybe Int
nBurn = Just 2000

nAutoTune :: Maybe Int
nAutoTune = Just 200

nIter :: Int
nIter = 20000

mhBench :: GenIO -> IO ()
mhBench g = do
  s <- mh nBurn nAutoTune nIter (mcmc 0 posterior moveCycle g)
  let t = map state . fromTrace $ trace s
      a = acceptance s
  putStrLn "Acceptance ratios:"
  putStrLn $ "Per move: " <> show (acceptanceRatios nIter a)
  putStrLn $ "Total: " <> show (acceptanceRatio nIter a)
  putStrLn "Mean and standard deviations:"
  putStrLn $ "True: " ++ show (trueMean, trueStdDev)
  putStrLn $ "Markov chain: " <> show (summarize t)
