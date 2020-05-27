{- |
Module      :  Poisson
Description :  Poisson regression model for airline fatalities
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May  8 13:30:51 2020.

See https://revbayes.github.io/tutorials/mcmc/poisson.html.

-}

module Poisson
  ( poissonBench
  ) where

import qualified Data.Vector.Unboxed as V
import Lens.Micro
import Numeric.Log hiding (sum)
import Statistics.Distribution hiding (mean, stdDev)
import Statistics.Distribution.Poisson
import Statistics.Sample
import System.Random.MWC

import Statistics.Mcmc

import Statistics.Mcmc.Acceptance
import Statistics.Mcmc.Item
import Statistics.Mcmc.Trace

type I = (Double, Double)

fatalities :: [Int]
fatalities = [24, 25, 31, 31, 22, 21, 26, 20, 16, 22]

normalizedYears :: [Double]
normalizedYears = map (subtract m) ys
  where
    ys = [1976.0 .. 1985.0]
    m = sum ys / fromIntegral (length ys)

f :: Int -> Double -> I -> Log Double
f ft yr (alpha, beta) = Exp $ logProbability (poisson l) (fromIntegral ft)
  where l = exp $ alpha + beta*yr

posterior :: I -> Log Double
posterior x = product [ f ft yr x | (ft, yr) <- zip fatalities normalizedYears ]

moveAlpha :: Move I
moveAlpha = slide "alpha" _1 0.0 0.2 True

moveBeta :: Move I
moveBeta = slide "beta" _2 0.0 0.2 True

moveCycle :: Cycle I
moveCycle = fromList [ (moveAlpha, 2)
                     , (moveBeta,  1) ]

initial :: I
initial = (0, 0)

summarize :: [I] -> ((Double, Double), (Double, Double))
summarize xs = ((mean as, stdDev as), (mean bs, stdDev bs))
  where as = V.fromList $ map fst xs
        bs = V.fromList $ map snd xs

nBurn :: Maybe Int
nBurn = Just 2000

nAutoTune :: Maybe Int
nAutoTune = Just 200

nIter :: Int
nIter = 10000

poissonBench :: GenIO -> IO ()
poissonBench g = do
  s <- mh nBurn nAutoTune nIter (mcmc initial posterior moveCycle g)
  let a = acceptance s
  putStrLn "Acceptance ratios:"
  putStrLn $ "Per move: " <> show (acceptanceRatios nIter a)
  putStrLn $ "Total: " <> show (acceptanceRatio nIter a)
  putStrLn "Mean and standard deviations:"
  let xs = map state . fromTrace $ trace s
      (ra, rb) = summarize xs
  putStrLn $ "Alpha: " <> show ra
  putStrLn $ "Beta: " <> show rb

