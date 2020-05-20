{- |
Module      :  Main
Description :  Benchmark with Criterion
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May  8 17:27:03 2020.

-}

module Main
  ( main
  ) where

import Criterion.Main
import System.Random.MWC

import Metropolis
import Poisson

main :: IO ()
main = do
  g <- create
  defaultMain
    [ bench "mh" $ nfIO (mhBench g)
    , bench "Poisson" $ nfIO (poissonBench g) ]
