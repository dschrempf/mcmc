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
import qualified Data.Vector.Storable as V
import Data.Vector.Storable (Vector)
import Lens.Micro
import Mcmc
import qualified Numeric.LinearAlgebra as L
import Numeric.LinearAlgebra (Herm, Matrix, R, (<#), (<.>))
import Numeric.Log
import Tree

-- TODO: Vectors. Ahh, we need to use a boxed data structure. For example, a
-- sequence. But then the likelihood function will be super slow.

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
pr xs | V.any (<= 0) xs = pzero
      | otherwise = Exp 0

-- Phylogenetic likelihood using a multivariate normal distribution. See
-- https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
--
-- The constant @k * log (2*pi)@ was left out on purpose.
--
-- lh meanVector invertedCovarianceMatrix logOfDeterminantOfCovarianceMatrix
lh :: Vector Double -> Matrix Double -> Double -> I -> Log Double
lh mu sigmaInv logSigmaDet xs = Exp $ ((-0.5) *) $ logSigmaDet + (dxs <# sigmaInv) <.> dxs
  where dxs = xs - mu

main :: IO ()
main = do
  trs <- manyNewick fn
  -- Get the posterior means and the posterior covariance matrix.
  let pm = getPosteriorMatrix trs
      (mu, sigma) = L.meanCov pm
      sigmaInv = L.inv $ L.unSym sigma
      logSigmaDet = log $ L.det $ L.unSym sigma
  -- Choose a bad starting state for our chain.
  let k = L.cols pm
      start = V.replicate k (1.0 :: Double)
  putStrLn "I am done."
