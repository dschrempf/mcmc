{- |
Module      :  Main
Description :  Benchmark Metropolis-Hastings algorithm
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 00:10:11 2020.

-}

module Main
  ( main
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
import Statistics.Mcmc.Move.Normal

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

posterior :: Double -> Log Double
posterior = Exp . logDensity (normalDistr trueMean trueStdDev)

moveCycle :: Cycle Double
moveCycle = Cycle [ moveNormalDouble "small" 0 0.1
                  , moveNormalDouble "medium" 0 1.0
                  , moveNormalDouble "large" 0 5.0
                  , moveNormalDouble "skewed" 1.0 4.0 ]

summarize :: [Double] -> (Double, Double)
summarize xs = (mean v, stdDev v)
  where v = V.fromList xs

-- TODO: This is more of a test, not a benchmark; use criterion.

main :: IO ()
main = do
  g <- create
  s <- mh 100000 (start 0 posterior moveCycle g)
  let t = map state . fromTrace $ trace s
  putStrLn "Acceptance ratios:"
  putStrLn $ "Per move: " <> show (map mvName $ fromCycle moveCycle) <> " " <> show (acceptanceRatios s)
  putStrLn $ "Total: " <> show (acceptanceRatio s)
  putStrLn "Mean and standard deviations:"
  putStrLn $ "True: " ++ show (trueMean, trueStdDev)
  putStrLn $ "Markov chain: " <> show (summarize t)
