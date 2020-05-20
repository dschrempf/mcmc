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
moveCycle = fromList [ (slideDouble "small"  0   0.1 False, 5)
                     , (slideDouble "medium" 0   1.0 False, 2)
                     , (slideDouble "large"  0   5.0 False, 2)
                     , (slideDouble "skewed" 1.0 4.0 False, 1)]

summarize :: [Double] -> (Double, Double)
summarize xs = (mean v, stdDev v)
  where v = V.fromList xs

n :: Int
n = 50000

mhBench :: GenIO -> IO ()
mhBench g = do
  s <- mh n (mcmc 0 posterior moveCycle g)
  let t = map state . fromTrace $ trace s
      a = acceptance s
  putStrLn "Acceptance ratios:"
  putStrLn $ "Per move: " <> show (acceptanceRatios n a)
  putStrLn $ "Total: " <> show (acceptanceRatio n a)
  putStrLn "Mean and standard deviations:"
  putStrLn $ "True: " ++ show (trueMean, trueStdDev)
  putStrLn $ "Markov chain: " <> show (summarize t)
