-- |
-- Module      :  Main
-- Copyright   :  (c) Dominik Schrempf 2020
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
import Normal
import Poisson
import System.Random.MWC

main :: IO ()
main = do
  g <- create
  defaultMain
    [ bgroup
        "Normal"
        [ bench "Slide" $ nfIO (normalSlide g),
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
        ]
    ]

-- benchmarking Normal
-- time                 43.83 ms   (43.57 ms .. 44.06 ms)
--                      1.000 R²   (1.000 R² .. 1.000 R²)
-- mean                 44.18 ms   (43.94 ms .. 44.55 ms)
-- std dev              615.6 μs   (311.7 μs .. 999.9 μs)

-- benchmarking NormalBactrian
-- time                 47.91 ms   (47.43 ms .. 48.58 ms)
--                      1.000 R²   (0.999 R² .. 1.000 R²)
-- mean                 48.14 ms   (47.91 ms .. 48.52 ms)
-- std dev              560.3 μs   (325.5 μs .. 734.8 μs)

-- benchmarking NormalCycle
-- time                 4.831 s    (4.796 s .. 4.855 s)
--                      1.000 R²   (1.000 R² .. 1.000 R²)
-- mean                 4.824 s    (4.813 s .. 4.829 s)
-- std dev              7.806 ms   (699.3 μs .. 9.704 ms)
-- variance introduced by outliers: 19% (moderately inflated)

-- benchmarking Poisson
-- time                 63.49 ms   (63.14 ms .. 64.09 ms)
--                      1.000 R²   (0.999 R² .. 1.000 R²)
-- mean                 63.50 ms   (63.37 ms .. 63.83 ms)
-- std dev              376.6 μs   (147.8 μs .. 702.4 μs)
