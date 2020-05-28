{- |
Module      :  Normal
Description :  Benchmark Metropolis-Hastings algorithm
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 00:10:11 2020.

-}

module Normal
  ( normalBench
  ) where

import qualified Data.Vector.Unboxed as V
import Numeric.Log as L
import Statistics.Distribution hiding (mean, stdDev)
import Statistics.Distribution.Normal
import Statistics.Sample
import System.Random.MWC

import Statistics.Mcmc

import Statistics.Mcmc.Item
import Statistics.Mcmc.Trace

trueMean :: Double
trueMean = 5

trueStdDev :: Double
trueStdDev = 4

posterior :: Double -> Log Double
posterior = Exp . logDensity (normalDistr trueMean trueStdDev)

moveCycle :: Cycle Double
moveCycle = fromList [ slideDouble "small"  5 0   0.1 True
                     , slideDouble "medium" 2 0   1.0 True
                     , slideDouble "large"  2 0   5.0 True
                     , slideDouble "skewed" 1 1.0 4.0 True ]

summarize :: [Double] -> (Double, Double)
summarize xs = (mean v, stdDev v)
  where v = V.fromList xs

nBurn :: Maybe Int
nBurn = Just 2000

nAutoTune :: Maybe Int
nAutoTune = Just 200

nIter :: Int
nIter = 20000

normalBench :: GenIO -> IO ()
normalBench g = do
  s <- mh nBurn nAutoTune nIter (mcmc 0 posterior moveCycle mempty g)
  let t = map state . fromTrace $ trace s
  putStrLn "Mean and standard deviations:"
  putStrLn $ "True: " ++ show (trueMean, trueStdDev)
  putStrLn $ "Markov chain: " <> show (summarize t)
