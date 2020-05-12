{- |
Module      :  Poisson
Description :  Poisson regression model for airline fatalities
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

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

import Statistics.Mcmc.Metropolis
import Statistics.Mcmc.Tools
import Statistics.Mcmc.Types
import Statistics.Mcmc.Move.Normal

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
moveAlpha = moveNormal _1 "alpha" 0.0 1.0

moveBeta :: Move I
moveBeta = moveNormal _2 "beta" 0.0 1.0

moveCycle :: Cycle I
moveCycle = Cycle [moveAlpha, moveBeta]

initial :: I
initial = (0, 0)

summarize :: [I] -> ((Double, Double), (Double, Double))
summarize xs = ((mean as, stdDev as), (mean bs, stdDev bs))
  where as = V.fromList $ map fst xs
        bs = V.fromList $ map snd xs

poissonBench :: GenIO -> IO ()
poissonBench g = do
  s <- mh 10000 (start initial posterior moveCycle g)
  putStrLn "Acceptance ratios:"
  putStrLn $ "Per move: " <> show (map mvName $ fromCycle moveCycle) <> " " <> show (acceptanceRatios s)
  putStrLn $ "Total: " <> show (acceptanceRatio s)
  putStrLn "Mean and standard deviations:"
  let xs = map state . fromTrace $ trace s
      -- ps = map (ln . logPosterior) . fromTrace $ trace s
      (ra, rb) = summarize xs
  -- print xs
  -- print ps
  putStrLn $ "Alpha: " <> show ra
  putStrLn $ "Beta: " <> show rb

