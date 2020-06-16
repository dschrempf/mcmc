{-# LANGUAGE OverloadedStrings #-}

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
  )
where

import qualified Data.Vector.Unboxed           as V
import           Lens.Micro
import           Numeric.Log             hiding ( sum )
import           Statistics.Distribution hiding ( mean
                                                , stdDev
                                                )
import           Statistics.Distribution.Poisson
import           Statistics.Sample
import           System.Random.MWC

import           Mcmc

import           Mcmc.Item
import           Mcmc.Trace

type I = (Double, Double)

fatalities :: [Int]
fatalities = [24, 25, 31, 31, 22, 21, 26, 20, 16, 22]

normalizedYears :: [Double]
normalizedYears = map (subtract m) ys
 where
  ys = [1976.0 .. 1985.0]
  m  = sum ys / fromIntegral (length ys)

f :: Int -> Double -> I -> Log Double
f ft yr (alpha, beta) = Exp $ logProbability (poisson l) (fromIntegral ft)
  where l = exp $ alpha + beta * yr

likelihood :: I -> Log Double
likelihood x =
  product [ f ft yr x | (ft, yr) <- zip fatalities normalizedYears ]

moveAlpha :: Move I
moveAlpha = slide "alpha" 2 _1 0.0 0.2 False

moveBeta :: Move I
moveBeta = slide "beta" 1 _2 0.0 0.2 False

moveCycle :: Cycle I
moveCycle = fromList [moveAlpha, moveBeta]

initial :: I
initial = (0, 0)

summarize :: [I] -> ((Double, Double), (Double, Double))
summarize xs = ((mean as, stdDev as), (mean bs, stdDev bs))
 where
  as = V.fromList $ map fst xs
  bs = V.fromList $ map snd xs

monAlpha :: MonitorParameter I
monAlpha = monitorRealFloat "alpha" _1

monBeta :: MonitorParameter I
monBeta = monitorRealFloat "beta" _2

monStd :: MonitorStdOut I
monStd = monitorStdOut [monAlpha, monBeta] 150

mon :: Monitor I
mon = Monitor monStd []

nBurn :: Maybe Int
nBurn = Just 2000

nAutoTune :: Maybe Int
nAutoTune = Just 200

nIter :: Int
nIter = 10000

poissonBench :: GenIO -> IO ()
poissonBench g = do
  let s = status (const 1) likelihood moveCycle mon initial g
  r <- mh nBurn nAutoTune nIter s
  putStrLn "Mean and standard deviations of:"
  let xs       = map state . fromTrace $ trace r
      (ra, rb) = summarize xs
  putStrLn $ "Alpha: " <> show ra
  putStrLn $ "Beta: " <> show rb
