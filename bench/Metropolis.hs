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

import Control.Monad
import Numeric.Log as L
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC

import Statistics.Mcmc.Metropolis
import Statistics.Mcmc.Tools
import Statistics.Mcmc.Types
import Statistics.Mcmc.Move.Normal

type I = Double

dDist :: NormalDistribution
dDist = normalDistr 5 0.3

pDist :: NormalDistribution
pDist = normalDistr 0 0.8

posterior :: [Double] -> I -> Log Double
posterior os x = L.sum [Exp $ logDensity pDist (x-o) | o <- os ]

moveCycle :: Cycle Double
moveCycle = Cycle [moveNormal "small" 0 0.1, moveNormal "medium" 0 0.2, moveNormal "large" 0 1.0]

-- TODO: This is more of a test, not a benchmark.

main :: IO ()
main = do
  g <- create
  os <- replicateM 60 (genContinuous dDist g)
  print os
  s <- mh 10000 (start 0 (posterior os) moveCycle g)
  print $ map state . fromTrace $ trace s
  print $ acceptanceRatios s
  print $ acceptanceRatio s
