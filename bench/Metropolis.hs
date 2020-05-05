{- |
Module      :  Main
Description :  Benchmark Metropolis-Hastings algorithm
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 00:10:11 2020.

TODO: Rename and create a main module (but is this necessary?).

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

type I = Double

post :: LogPosterior I
post (S x) = Exp $ log $ density (normalDistr 0 0.3) x

start :: Item I
start = Item (S 0) (post (S 0))

mvDist :: NormalDistribution
mvDist = normalDistr 0 0.1

mvS :: S I -> GenIO -> IO (S I)
mvS (S x) g = do
  d <- genContinuous mvDist g
  return $ S (x + d)

move :: Move I
move = Move mvS (\(S x) (S y) -> Exp $ log $ density mvDist (y-x))

status :: GenIO -> Status I
status = Status start post (Cycle [move]) (Trace [start])

main :: IO ()
main = do
  g <- createSystemRandom
  s <- mh 500 (status g)
  let (Trace is) = mcmcTrace s
  print $ map (unpack . lnState) is
