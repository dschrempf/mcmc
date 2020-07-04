{-# LANGUAGE BangPatterns #-}

-- |
-- Module      :  Main
-- Description :  Approximate phylogenetic likelihood with multivariate normal distribution
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri Jul  3 21:27:38 2020.
--
-- This example is a little more involved and includes estimation of a larger
-- phylogeny by fitting the posterior with a multivariate normal distribution.
--
-- The trees are read from a data file which is given relative to the directory
-- in which this module resides. Hence, the compiled binary has to be executed
-- from this directory.
module Main (main) where

import Algebra.Graph.Labelled.AdjacencyMap
import Control.Monad
import qualified Data.Vector.Storable as V
import Data.Vector.Storable (Vector)
import Lens.Micro.Platform
import Mcmc
import qualified Numeric.LinearAlgebra as L
import Numeric.LinearAlgebra ((<#), (<.>), Matrix, R)
import Numeric.Log
import System.Random.MWC
import Tree

-- TODO: Move on vectors. This is really hard.

-- TODO: Monitors.

-- We condense the branch lengths into a vector.
type I = Vector R

fn :: FilePath
fn = "data/plants_1.treelist.gz"

getEdges :: T Double a -> I
getEdges = L.fromList . (map (fromD . (^. _1)) . edgeList)

getPosteriorMatrix :: [T Double a] -> Matrix R
getPosteriorMatrix = L.fromRows . map getEdges

-- Uniform prior. Ensuring positive branch lengths. If this is too slow, the
-- positiveness of branches has to be ensured by the moves.
pr :: I -> Log Double
pr xs
  | V.any (<= 0) xs = pzero
  | otherwise = Exp 0

-- Phylogenetic likelihood using a multivariate normal distribution. See
-- https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
--
-- The constant @k * log (2*pi)@ was left out on purpose.
--
-- lh meanVector invertedCovarianceMatrix logOfDeterminantOfCovarianceMatrix
lh :: Vector Double -> Matrix Double -> Double -> I -> Log Double
lh mu sigmaInv logSigmaDet xs = Exp $ ((-0.5) *) $ logSigmaDet + (dxs <# sigmaInv) <.> dxs
  where
    dxs = xs - mu

-- Slide branch with given index.
slideBranch :: Int -> Move I
slideBranch i = slideSymmetric n 1 (singular $ ix i) 0.1 True
  where
    n = "Slide branch " <> show i

moveCycle :: I -> Cycle I
moveCycle v = fromList [slideBranch i | i <- [0 .. k]]
  where
    k = V.length v - 1

-- Branch length monitors.
branchMons :: I -> [MonitorParameter I]
branchMons v = [monitorRealFloat (n i) (singular $ ix i) | i <- [0 .. k]]
  where
    n i = "Branch " <> show i
    k = V.length v - 1

mon :: I -> Monitor I
mon v = Monitor (monitorStdOut (take 3 bs) 50) [monitorFile "Branches" bs 10] []
  where bs = branchMons v

-- Number of burn in iterations.
nBurnIn :: Maybe Int
nBurnIn = Just 2000

-- Auto tuning period.
nAutoTune :: Maybe Int
nAutoTune = Just 100

-- Number of Metropolis-Hasting iterations after burn in.
nIterations :: Int
nIterations = 10000

main :: IO ()
main = do
  g <- create
  putStrLn "Read trees."
  !trs <- nNewick 300 fn
  putStrLn "Get the posterior means and the posterior covariance matrix."
  let !pm = getPosteriorMatrix trs
      !(mu, sigma) = L.meanCov pm
      !sigmaInv = L.inv $ L.unSym sigma
      !logSigmaDet = log $ L.det $ L.unSym sigma
  print $ V.length mu
  print $ L.rows $ L.unSym sigma
  -- TODO: Infinity?
  print logSigmaDet
  putStrLn "Choose a bad starting state for our chain."
  let k = L.cols pm
      start = V.replicate k (1.0 :: Double)
  -- Status of the chain.
  let s = debug $ status "ApproximatePhylogeneticLikelihoodMultivariate" pr (lh mu sigmaInv logSigmaDet) (moveCycle start) (mon start) start nBurnIn nAutoTune nIterations g
  void $ mh s
  putStrLn "I am done."
