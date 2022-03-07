-- |
-- Module      :  Main
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May  8 17:27:03 2020.
module Main
  ( main,
  )
where

import Criterion.Main
import Data.List
import Mcmc.Internal.Gamma
import Normal
import Numeric.SpecFunctions
import Poisson
import System.Random.MWC

gammaBenchG :: RealFloat a => (a -> a) -> [a] -> a
gammaBenchG f = foldl' (\acc x -> acc + (f x)) 0
{-# SPECIALIZE gammaBenchG :: (Double -> Double) -> [Double] -> Double #-}

gammaVals :: [Double]
gammaVals = [0, 0.01 .. 10000]

main :: IO ()
main = do
  g <- create
  defaultMain
    [ bgroup
        "Normal"
        [ bench "Slide" $ nfIO (normalSlideBench g),
          bench "Bactrian" $ nfIO (normalBactrianBench g),
          bench "LargeCycle" $ nfIO (normalLargeCycleBench g)
        ],
      bench "Poisson" $ nfIO (poissonBench g),
      bgroup
        "MC3"
        [ bench "MC3 2" $ nfIO (normalMC3 g 2),
          bench "MC3 3" $ nfIO (normalMC3 g 3),
          bench "MC3 4" $ nfIO (normalMC3 g 4),
          bench "MC3 5" $ nfIO (normalMC3 g 5),
          bench "MC3 10" $ nfIO (normalMC3 g 10)
        ],
      bgroup
        "probability"
        [ bench "gamma function general" $
            nf (gammaBenchG logGammaG) gammaVals,
          bench "gamma function specialized" $
            nf (gammaBenchG logGamma) gammaVals
        ]
    ]

-- benchmarking Normal/Slide
-- time                 42.31 ms   (41.88 ms .. 42.60 ms)
--                      1.000 R²   (0.999 R² .. 1.000 R²)
-- mean                 42.75 ms   (42.52 ms .. 43.29 ms)
-- std dev              661.4 μs   (347.4 μs .. 1.074 ms)

-- benchmarking Normal/Bactrian
-- time                 45.51 ms   (45.30 ms .. 45.92 ms)
--                      1.000 R²   (0.999 R² .. 1.000 R²)
-- mean                 45.41 ms   (45.31 ms .. 45.61 ms)
-- std dev              276.0 μs   (141.5 μs .. 460.8 μs)

-- benchmarking Normal/LargeCycle
-- time                 68.82 ms   (67.18 ms .. 70.81 ms)
--                      0.999 R²   (0.997 R² .. 1.000 R²)
-- mean                 67.68 ms   (67.26 ms .. 68.59 ms)
-- std dev              1.074 ms   (618.2 μs .. 1.602 ms)

-- benchmarking Poisson
-- time                 72.94 ms   (63.07 ms .. 87.73 ms)
--                      0.953 R²   (0.920 R² .. 1.000 R²)
-- mean                 64.76 ms   (62.84 ms .. 71.62 ms)
-- std dev              5.785 ms   (783.9 μs .. 10.08 ms)
-- variance introduced by outliers: 26% (moderately inflated)

-- benchmarking MC3/MC3 2
-- time                 13.08 ms   (12.73 ms .. 13.44 ms)
--                      0.993 R²   (0.986 R² .. 0.997 R²)
-- mean                 13.41 ms   (13.16 ms .. 13.72 ms)
-- std dev              682.7 μs   (520.3 μs .. 874.0 μs)
-- variance introduced by outliers: 22% (moderately inflated)

-- benchmarking MC3/MC3 3
-- time                 19.19 ms   (18.86 ms .. 19.59 ms)
--                      0.998 R²   (0.996 R² .. 1.000 R²)
-- mean                 19.28 ms   (19.11 ms .. 19.51 ms)
-- std dev              454.1 μs   (339.1 μs .. 608.4 μs)

-- benchmarking MC3/MC3 4
-- time                 25.01 ms   (24.21 ms .. 25.66 ms)
--                      0.997 R²   (0.996 R² .. 0.999 R²)
-- mean                 24.21 ms   (23.99 ms .. 24.55 ms)
-- std dev              606.1 μs   (414.7 μs .. 738.5 μs)

-- benchmarking MC3/MC3 5
-- time                 28.39 ms   (26.99 ms .. 29.34 ms)
--                      0.995 R²   (0.990 R² .. 0.999 R²)
-- mean                 31.13 ms   (30.09 ms .. 33.49 ms)
-- std dev              3.009 ms   (984.0 μs .. 4.319 ms)
-- variance introduced by outliers: 40% (moderately inflated)

-- benchmarking MC3/MC3 10
-- time                 57.25 ms   (56.98 ms .. 57.61 ms)
--                      1.000 R²   (1.000 R² .. 1.000 R²)
-- mean                 57.46 ms   (57.34 ms .. 57.56 ms)
-- std dev              192.4 μs   (140.1 μs .. 284.0 μs)
