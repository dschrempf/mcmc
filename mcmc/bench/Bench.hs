-- |
-- Module      :  Main
-- Copyright   :  2021 Dominik Schrempf
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
import Mcmc.Internal.SpecFunctions
import Normal
import Numeric.SpecFunctions
import Poisson
import System.Random.Stateful

gammaBenchG :: RealFloat a => (a -> a) -> [a] -> a
gammaBenchG f = foldl' (\acc x -> acc + f x) 0
{-# SPECIALIZE gammaBenchG :: (Double -> Double) -> [Double] -> Double #-}

gammaVals :: [Double]
gammaVals = [0, 0.01 .. 10000]

main :: IO ()
main = do
  let g = mkStdGen 0
  defaultMain
    [ bgroup
        "Normal"
        [ bench "Slide" $ nfIO (normalSlideBench g),
          bench "Bactrian" $ nfIO (normalBactrianBench g),
          bench "LargeCycle" $ nfIO (normalLargeCycleBench g)
        ],
      bgroup
        "Poisson"
        [ bench "Random walk" $ nfIO (poissonBench g),
          bench "Hamiltonian" $ nfIO (poissonHamiltonianBench g)
        ],
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
-- time                 24.12 ms   (23.52 ms .. 24.76 ms)
--                      0.997 R²   (0.995 R² .. 0.999 R²)
-- mean                 23.70 ms   (23.44 ms .. 24.02 ms)
-- std dev              658.2 μs   (448.5 μs .. 837.2 μs)

-- benchmarking Normal/Bactrian
-- time                 26.65 ms   (26.08 ms .. 27.26 ms)
--                      0.998 R²   (0.996 R² .. 1.000 R²)
-- mean                 26.53 ms   (26.32 ms .. 26.91 ms)
-- std dev              626.9 μs   (322.3 μs .. 1.040 ms)

-- benchmarking Normal/LargeCycle
-- time                 52.48 ms   (51.51 ms .. 53.83 ms)
--                      0.999 R²   (0.998 R² .. 1.000 R²)
-- mean                 52.54 ms   (51.62 ms .. 54.90 ms)
-- std dev              2.765 ms   (758.0 μs .. 4.788 ms)
-- variance introduced by outliers: 14% (moderately inflated)

-- benchmarking Poisson/Random walk
-- time                 86.57 ms   (86.06 ms .. 86.97 ms)
--                      1.000 R²   (1.000 R² .. 1.000 R²)
-- mean                 86.57 ms   (86.20 ms .. 87.14 ms)
-- std dev              775.9 μs   (415.4 μs .. 1.205 ms)

-- benchmarking Poisson/Hamiltonian
-- time                 2.890 s    (2.313 s .. 3.359 s)
--                      0.993 R²   (0.992 R² .. 1.000 R²)
-- mean                 3.082 s    (2.960 s .. 3.153 s)
-- std dev              120.4 ms   (39.93 ms .. 164.0 ms)
-- variance introduced by outliers: 19% (moderately inflated)

-- benchmarking MC3/MC3 2
-- time                 7.771 ms   (7.573 ms .. 8.008 ms)
--                      0.996 R²   (0.993 R² .. 0.999 R²)
-- mean                 7.810 ms   (7.722 ms .. 7.974 ms)
-- std dev              328.3 μs   (223.0 μs .. 517.6 μs)
-- variance introduced by outliers: 19% (moderately inflated)

-- benchmarking MC3/MC3 3
-- time                 10.69 ms   (10.57 ms .. 10.87 ms)
--                      0.997 R²   (0.993 R² .. 1.000 R²)
-- mean                 10.73 ms   (10.65 ms .. 10.90 ms)
-- std dev              308.1 μs   (132.2 μs .. 528.9 μs)

-- benchmarking MC3/MC3 4
-- time                 13.81 ms   (13.66 ms .. 14.08 ms)
--                      0.999 R²   (0.997 R² .. 1.000 R²)
-- mean                 13.73 ms   (13.67 ms .. 13.85 ms)
-- std dev              213.2 μs   (95.38 μs .. 334.4 μs)

-- benchmarking MC3/MC3 5
-- time                 16.80 ms   (16.65 ms .. 16.93 ms)
--                      0.999 R²   (0.997 R² .. 1.000 R²)
-- mean                 16.99 ms   (16.87 ms .. 17.23 ms)
-- std dev              376.6 μs   (162.3 μs .. 589.3 μs)

-- benchmarking MC3/MC3 10
-- time                 32.59 ms   (32.10 ms .. 33.20 ms)
--                      0.999 R²   (0.998 R² .. 1.000 R²)
-- mean                 32.21 ms   (32.07 ms .. 32.46 ms)
-- std dev              390.9 μs   (219.6 μs .. 640.0 μs)

-- benchmarking probability/gamma function general
-- time                 60.73 ms   (60.37 ms .. 61.26 ms)
--                      1.000 R²   (1.000 R² .. 1.000 R²)
-- mean                 60.70 ms   (60.49 ms .. 60.91 ms)
-- std dev              385.5 μs   (287.3 μs .. 513.5 μs)

-- benchmarking probability/gamma function specialized
-- time                 59.00 ms   (58.74 ms .. 59.32 ms)
--                      1.000 R²   (1.000 R² .. 1.000 R²)
-- mean                 58.83 ms   (58.76 ms .. 58.97 ms)
-- std dev              168.3 μs   (105.7 μs .. 245.7 μs)
