{-# LANGUAGE OverloadedStrings #-}

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
-- This example is more involved and includes estimation of a larger phylogeny
-- by fitting the posterior with a multivariate normal distribution.
--
-- The trees are read from a file which is given relative to the @mcmc@ git
-- repository base directory. Hence, the compiled binary has to be executed from
-- the base directory.
--
-- The module hierarchy is organized as follows:
--
-- - Main: Functions to prepare the data, run and continue the
--   Metropolis-Hasting sampler, and to inspect the application.
--
-- - Definitions: The state space, prior distribution, and the likelihood
--   function of the sampler. Also includes the proposals and the monitor.
--
-- - Calibration and Constrain: Calibrations on node ages and node order
--   constraints.
--
-- - Tools: Miscellaneous tools to prepare the data.
--
-- If you try to understand what is going on, your starting point should be
-- 'Definitions'.
module Main
  ( main,
  )
where

import Control.Lens
import Control.Monad
import Criterion
import Data.Aeson
import Data.Bifunctor
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List
import Data.Maybe
import qualified Data.Set as S
import qualified Data.Vector.Storable as V
import qualified ELynx.Topology.Rooted as T
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import System.Environment
import System.Random.MWC hiding (uniform)

-- Disable the syntax formatter Ormolu to highlight relevant module imports.
{- ORMOLU_DISABLE -}
-- The ELynx library includes functions to work on trees.
import ELynx.Tree

-- The Mcmc library includes the Metropolis-Hastings sampler.
import Mcmc
import Mcmc.Tree

-- Local modules (see comment above).
import Calibrations
import Constraints
import Definitions
import Tools
{- ORMOLU_ENABLE -}

fnMeanTree :: FilePath
fnMeanTree = bnAnalysis ++ ".meantree"

-- The rooted tree with posterior mean branch lengths will be stored in a file
-- with this name.
getMeanTree :: IO (Tree Double BS.ByteString)
getMeanTree = oneTree fnMeanTree

fnData :: FilePath
fnData = bnAnalysis ++ ".data"

-- Get the posterior branch length means, the inverted covariance matrix, and
-- the determinant of the covariance matrix.
getData :: IO (V.Vector Double, L.Matrix Double, Double)
getData = do
  (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' fnData
  let sigmaInv = L.fromRows sigmaInvRows
  return (mu, sigmaInv, logSigmaDet)

-- Get the posterior matrix of branch lengths of rooted trees. Merge the two
-- branch lengths leading to the root.
getPosteriorMatrixRooted :: [Tree Double a] -> L.Matrix Double
getPosteriorMatrixRooted = L.fromRows . map (sumFirstTwo . getBranches)

-- Get the posterior matrix of branch lengths of unrooted trees (trees with
-- multifurcating root nodes). Before midpoint rooting, the mean branch lengths
-- of the unrooted trees have to be determined.
getPosteriorMatrix :: [Tree Double a] -> L.Matrix Double
getPosteriorMatrix = L.fromRows . map (V.fromList . branches)

-- Only use this if absolutely necessary...
beautifyVariance :: Double -> Double -> Double
beautifyVariance _ x
  | x < 0 = error "beautifyVariance: Variance is negative."
  | x < eps = eps
  | otherwise = x
  where
    eps = 1e-4

-- Analyze the covariance matrix and change problematic values. This step is
-- awful but necessary when there are not enough samples from the posterior.
beautifyCovarianceMatrix :: L.Matrix Double -> L.Matrix Double
beautifyCovarianceMatrix m =
  L.accum m beautifyVariance [((i, i), 0) | i <- [0 .. nRows - 1]]
  where
    nRows = L.rows m

-- Read trees and extract branch lengths.
prepare :: IO ()
prepare = do
  putStrLn "Read trees."
  trsAll <- someTrees Standard fnInTrees
  let nTrees = length trsAll
      nBurnInTrees = nTrees `div` 10
      trs = drop nBurnInTrees trsAll
  putStrLn $ show nTrees ++ " read; skip a burn in of " ++ show nBurnInTrees ++ " trees."

  putStrLn "Check if trees have the same topology."
  let l = length $ nub $ map T.fromLabeledTree trs
  unless (l == 1) (error "Trees have different topologies.")

  putStrLn "Calculate mean branch lengths."
  let pm = getPosteriorMatrix trs
      (means, _) = L.meanCov pm
  putStrLn "The mean branch lengths are:"
  print means
  let meanTree =
        fromMaybe (error "Could not label tree with mean branch lengths") $
          setBranches (map Length $ V.toList means) (head trs)
      lvs = leaves meanTree
      trOutgroup = either error id $ outgroup (S.singleton $ head lvs) "root" meanTree
      tr = either error id $ midpoint trOutgroup
  putStrLn "The tree with mean branch lengths rooted at the midpoint:"
  print $ toNewick $ measurableToPhyloTree tr
  putStrLn $ "Save the mean tree to " <> fnMeanTree <> "."
  BL.writeFile fnMeanTree (toNewick $ measurableToPhyloTree tr)

  putStrLn "Root the trees at the midpoint of the mean tree."
  let outgroups = fst $ fromBipartition $ either error id $ bipartition tr
      trsRooted = map (first fromLength . either error id . outgroup outgroups "root" . first Length) trs
  putStrLn "Get the posterior means and the posterior covariance matrix."
  let pmR = getPosteriorMatrixRooted trsRooted
      (mu, sigmaBare) = second L.unSym $ L.meanCov pmR
  putStrLn $ "Minimum value of absolute values of covariance matrix: " ++ show (L.minElement $ L.cmap abs sigmaBare)
  putStrLn $ "Maximum value of absolute values of covariance matrix: " ++ show (L.maxElement $ L.cmap abs sigmaBare)
  let variancesBare = L.takeDiag sigmaBare
  putStrLn "The variances are: "
  print variancesBare
  putStrLn $ "Minimum variance: " ++ show (L.minElement variancesBare)
  putStrLn $ "Maximum variance: " ++ show (L.maxElement variancesBare)

  putStrLn "Beautify covariance matrix. Ouch!"
  let sigma = beautifyCovarianceMatrix sigmaBare
      variances = L.takeDiag sigma
  putStrLn $ "Minimum variance: " ++ show (L.minElement variances)
  putStrLn $ "Maximum variance: " ++ show (L.maxElement variances)

  putStrLn "Prepare the covariance matrix for the likelihood calculation."
  let (sigmaInv, (logSigmaDet, _)) = L.invlndet sigma
  putStrLn $ "The logarithm of the determinant of the covariance matrix is: " ++ show logSigmaDet

  putStrLn $ "Save the posterior means and covariances to " <> fnData <> "."
  encodeFile fnData (mu, L.toRows sigmaInv, logSigmaDet)

-- Run the Metropolis-Hastings sampler.
runMetropolisHastings :: IO ()
runMetropolisHastings = do
  -- Read the mean tree and the posterior means and covariances.
  meanTree <- getMeanTree
  (mu, sigmaInv, logSigmaDet) <- getData
  putStrLn "The posterior means of the branch lengths are:"
  print mu
  -- Initialize a starting state using the mean tree.
  let start = initWith $ identify meanTree
      cb = getCalibrations meanTree
      cs = getConstraints meanTree
      pr' = priorDistribution cb cs
      lh' = likelihoodFunction mu sigmaInv logSigmaDet
      ccl' = proposals $ identify meanTree
  -- Create a seed value for the random number generator. Actually, the
  -- 'create' function is deterministic, but useful during development. For
  -- real analyses, use 'createSystemRandom'.
  g <- create
  -- Construct the status of the Markov chain.
  let s =
        force . cleanWith cleaner . saveWith 1 . debug $
          -- Have a look at the 'status' function to understand the
          -- different parameters.
          status
            bnAnalysis
            pr'
            lh'
            ccl'
            (monitor cb cs)
            start
            nBurnIn
            nAutoTune
            nIterations
            g
  -- Run the Markov chain.
  void $ mh s

continueMetropolisHastings :: Int -> IO ()
continueMetropolisHastings n = do
  meanTree <- getMeanTree
  (mu, sigmaInv, logSigmaDet) <- getData
  let cb = getCalibrations meanTree
      cs = getConstraints meanTree
  -- Load the MCMC status.
  s <-
    loadStatus
      (priorDistribution cb cs)
      (likelihoodFunction mu sigmaInv logSigmaDet)
      (Just cleaner)
      (proposals $ identify meanTree)
      (monitor cb cs)
      (bnAnalysis ++ ".mcmc")
  void $ mhContinue n s

convert :: IO ()
convert = do
  meanTree <- getMeanTree
  let tbl = zip (map BS.unpack $ leaves meanTree) (map show $ leaves $ identify meanTree)
  putStr $ unlines $ map show tbl

-- Benchmark different functions used by the MCMC sampler.
runBench :: IO ()
runBench = do
  tr <- getMeanTree
  let (pth, _) =
        fromMaybe (error "Gn_montanu not found.") $
          ifind (\_ n -> n == "Gn_montanu") tr
  putStrLn $ "The path to \"Gn_montanu\" is: " <> show pth
  let bf1 =
        toTree . insertLabel "Bla"
          . fromMaybe (error "Path does not lead to a leaf.")
          . goPath pth
          . fromTree
  putStrLn $ "Change a leaf: " <> show (bf1 tr) <> "."
  putStrLn "Benchmark change a leaf."
  benchmark $ nf bf1 tr
  let bf2 =
        label . current
          . fromMaybe (error "Path does not lead to a leaf.")
          . goPath pth
          . fromTree
  putStrLn $ "Leaf to get: " <> show (bf2 tr) <> "."
  putStrLn "Benchmark get a leaf."
  benchmark $ nf bf2 tr
  putStrLn "Benchmark calculation of prior."
  let i = initWith $ identify tr
      pr' = priorDistribution (getCalibrations tr) (getConstraints tr)
  benchmark $ nf pr' i
  (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' fnData
  let sigmaInv = L.fromRows sigmaInvRows
      lh' = likelihoodFunction mu sigmaInv logSigmaDet
  putStrLn "Benchmark calculation of likelihood."
  benchmark $ nf lh' i

  putStrLn "Benchmark identify."
  benchmark $ nf identify tr

-- Inspect different objects; useful for debugging.
inspect :: IO ()
inspect = do
  tr <- getMeanTree
  let lvs = leaves tr
  putStrLn $ "The mean tree has " <> show (length lvs) <> " leaves."

  let i = initWith $ identify tr
      pr' = priorDistribution (getCalibrations tr) (getConstraints tr)
  putStrLn $ "Test if time tree is ultrametric: " <> show (ultrametric $ _timeTree i)
  putStrLn $ "Initial prior: " <> show (pr' i) <> "."
  putStrLn "Load posterior means and covariances."
  (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' fnData
  let sigmaInv = L.fromRows sigmaInvRows
      lh' = likelihoodFunction mu sigmaInv logSigmaDet
  putStrLn $ "Initial log-likelihood: " <> show (ln $ lh' i) <> "."

  putStrLn "Paths to calibration nodes."
  putStrLn "Slight difference: 179, 183 younger than 245 247."
  putStrLn $ "Mrca of 179 and 183: " ++ show (mrca [179, 183] $ identify tr)
  putStrLn $ "Mrca of 245 and 247: " ++ show (mrca [245, 247] $ identify tr)
  putStrLn "More extreme difference: 34, 38 younger than 62, 64."
  putStrLn $ "Mrca of 34 and 38: " ++ show (mrca [34, 38] $ identify tr)
  putStrLn $ "Mrca of 62 and 64: " ++ show (mrca [62, 64] $ identify tr)

main :: IO ()
main = do
  -- Get arguments.
  as <- getArgs
  case as of
    -- Read in all trees, calculate posterior means and covariances of the
    -- branch lengths, and find the midpoint root of the mean tree.
    ["prepare"] -> do
      prepare
    -- Run the Metropolis-Hastings sampler.
    ["run"] -> do
      runMetropolisHastings
    -- Continue sampling.
    ["continue", n] -> do
      continueMetropolisHastings (read n)
    -- Print conversion of leaves.
    ["convert"] -> do
      convert
    -- Benchmark different functions used by the MCMC sampler.
    ["bench"] -> do
      runBench
    -- Inspect different objects; useful for debugging.
    ["inspect"] -> do
      inspect
    -- Print usage instructions if none of the previous commands was entered.
    _ -> putStrLn "Use one command of: [prepare|run|continue N|convert|bench|inspect]."

-- Post scriptum:
--
-- Benchmarks with criterion indicate the following:
--
-- 1. Algebraic graphs are slow, but the adjacency map with 'Int' node labels and branch labels is reasonably fast.
-- 2. Lenses are fast, but don't allow access children, or the parent.
-- 3. Zippers are very slow.
--
-- @
-- AdjacencyIntMap; postSet 150
-- benchmarking...
-- time                 41.25 ns   (41.09 ns .. 41.58 ns)
--                      1.000 R²   (0.999 R² .. 1.000 R²)
-- mean                 41.21 ns   (41.13 ns .. 41.45 ns)
-- std dev              467.6 ps   (126.7 ps .. 873.2 ps)
-- variance introduced by outliers: 12% (moderately inflated)
--
-- AdjacencyIntMap; preSet 150
-- benchmarking...
-- time                 5.335 μs   (5.308 μs .. 5.372 μs)
--                      0.999 R²   (0.998 R² .. 1.000 R²)
-- mean                 5.335 μs   (5.310 μs .. 5.450 μs)
-- std dev              146.3 ns   (31.19 ns .. 323.8 ns)
-- variance introduced by outliers: 33% (moderately inflated)
--
-- AdjacencyIntMap; replaceEdge 1 12
-- benchmarking...
-- time                 5.693 μs   (5.617 μs .. 5.819 μs)
--                      0.994 R²   (0.985 R² .. 0.999 R²)
-- mean                 5.816 μs   (5.671 μs .. 6.089 μs)
-- std dev              620.6 ns   (261.2 ns .. 1.094 μs)
-- variance introduced by outliers: 89% (severely inflated)
--
-- Control.Zipper; preview node 150
-- benchmarking...
-- time                 2.513 μs   (2.351 μs .. 2.678 μs)
--                      0.978 R²   (0.972 R² .. 0.993 R²)
-- mean                 2.432 μs   (2.345 μs .. 2.546 μs)
-- std dev              336.4 ns   (248.6 ns .. 467.2 ns)
-- variance introduced by outliers: 93% (severely inflated)
--
-- Control.Zipper; set node 150
-- benchmarking...
-- time                 8.838 μs   (8.768 μs .. 8.937 μs)
--                      0.999 R²   (0.998 R² .. 1.000 R²)
-- mean                 8.815 μs   (8.770 μs .. 8.895 μs)
-- std dev              205.4 ns   (116.7 ns .. 306.8 ns)
-- variance introduced by outliers: 25% (moderately inflated)
--
-- Specialized zipper; read a leaf node label
-- benchmarking...
-- time                 348.7 ns   (320.3 ns .. 386.3 ns)
--                      0.947 R²   (0.933 R² .. 0.974 R²)
-- mean                 351.8 ns   (330.0 ns .. 373.6 ns)
-- std dev              75.84 ns   (58.43 ns .. 88.84 ns)
-- variance introduced by outliers: 98% (severely inflated)
--
-- Specialized zipper; change a leaf node (this is not possible with
-- Control.Zipper).
-- benchmarking...
-- time                 5.389 μs   (5.382 μs .. 5.401 μs)
--                      0.999 R²   (0.999 R² .. 1.000 R²)
-- mean                 5.466 μs   (5.386 μs .. 5.686 μs)
-- std dev              415.0 ns   (15.48 ns .. 773.6 ns)
-- variance introduced by outliers: 79% (severely inflated)
--
--
-- Given the last benchmark, this would mean that for 500 changes per iteration
-- and 100k iterations, we spend 300 seconds traversing and changing the tree
-- without calculating anything else. I guess that's OK.
-- @
