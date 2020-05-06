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

import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC

import Statistics.Mcmc.Types
import Statistics.Mcmc.Metropolis
import Statistics.Mcmc.Move.Normal

type I = Double

posterior :: I -> Log Double
posterior x = Exp $ log $ density (normalDistr 0 0.3) x

moveCycle :: Cycle Double
moveCycle = Cycle [moveNormal 0 0.1, moveNormal 0 0.2]

main :: IO ()
main = do
  g <- create
  s <- mh 5 (start 0 posterior moveCycle g)
  let (Trace is) = trace s
  print $ map state is
