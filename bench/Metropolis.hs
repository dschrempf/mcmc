{- |
Module      :  Metropolis
Description :  Benchmark Metropolis-Hastings algorithm
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

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

import Statistics.Mcmc.Metropolis
import Statistics.Mcmc.Tools
import Statistics.Mcmc.Types
import Statistics.Mcmc.Moves

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

posterior :: Double -> Log Double
posterior = Exp . logDensity (normalDistr trueMean trueStdDev)

moveCycle :: Cycle Double
moveCycle = Cycle [ slideDouble "small" 0 0.1
                  , slideDouble "medium" 0 1.0
                  , slideDouble "large" 0 5.0
                  , slideDouble "skewed" 1.0 4.0 ]

summarize :: [Double] -> (Double, Double)
summarize xs = (mean v, stdDev v)
  where v = V.fromList xs

mhBench :: GenIO -> IO ()
mhBench g = do
  s <- mh 50000 (start 0 posterior moveCycle g)
  let t = map state . fromTrace $ trace s
  putStrLn "Acceptance ratios:"
  putStrLn $ "Per move: " <> show (map mvName $ fromCycle moveCycle) <> " " <> show (acceptanceRatios s)
  putStrLn $ "Total: " <> show (acceptanceRatio s)
  putStrLn "Mean and standard deviations:"
  putStrLn $ "True: " ++ show (trueMean, trueStdDev)
  putStrLn $ "Markov chain: " <> show (summarize t)
