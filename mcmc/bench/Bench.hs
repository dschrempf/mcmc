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
    [ bench "Normal" $ nfIO (normalBench g),
      bench "NormalBactrian" $ nfIO (normalBactrianBench g),
      bench "Poisson" $ nfIO (poissonBench g)
    ]
