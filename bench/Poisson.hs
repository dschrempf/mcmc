{- |
Module      :  Main
Description :  Poisson regression model for airline fatalities
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May  8 13:30:51 2020.

See https://revbayes.github.io/tutorials/mcmc/poisson.html.

-}

module Main
  ( main
  ) where

import qualified Data.Vector.Unboxed as V
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

years :: [Int]
years = [1976 .. 1985]

normalizedYears :: [Double]
normalizedYears = map (subtract m . fromIntegral) years
  where
    m :: Double
    m = fromIntegral (sum years) / fromIntegral (length years)

pDist :: Int -> Double -> I -> Log Double
pDist ft yr (alpha, beta) = Exp $ logProbability (poisson l) (fromIntegral ft)
  where l = exp $ alpha + beta*yr

posterior :: I -> Log Double
posterior x = product [ pDist ft yr x | (ft, yr) <- zip fatalities normalizedYears ]

moveAlpha :: Move I
moveAlpha = moveNormal fst (\a (_, b) -> (a, b)) "alpha" 0.0 1.0

moveBeta :: Move I
moveBeta = moveNormal snd (\b (a, _) -> (a, b)) "beta" 0.0 1.0

moveCycle :: Cycle I
moveCycle = Cycle [moveAlpha, moveBeta]

initial :: I
initial = (0, 0)

summarize :: [I] -> ((Double, Double), (Double, Double))
summarize xs = ((mean as, stdDev as), (mean bs, stdDev bs))
  where as = V.fromList $ map fst xs
        bs = V.fromList $ map snd xs

-- TODO: This is more of a test, not a benchmark; use criterion.

main :: IO ()
main = do
  g <- create
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

