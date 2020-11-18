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

import Control.Monad
import Data.Aeson
import Data.Bifunctor
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List
import Data.Maybe
import qualified Data.Set as S
import qualified Data.Vector.Storable as V
import qualified Numeric.LinearAlgebra as L
import System.Environment
import System.Random.MWC hiding (uniform)

-- Disable the syntax formatter Ormolu to highlight relevant module imports.
{- ORMOLU_DISABLE -}
-- The ELynx library includes functions to work on trees.
import qualified ELynx.Topology as T
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
getMeanTree :: IO (Tree Length Name)
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
getPosteriorMatrixRooted :: [Tree Length a] -> L.Matrix Double
getPosteriorMatrixRooted = L.fromRows . map (sumFirstTwo . getBranches)

-- Get the posterior matrix of branch lengths of unrooted trees (trees with
-- multifurcating root nodes). Before midpoint rooting, the mean branch lengths
-- of the unrooted trees have to be determined.
getPosteriorMatrix :: [Tree Length a] -> L.Matrix Double
getPosteriorMatrix = L.fromRows . map (V.fromList . map fromLength . branches)

-- -- Only use this if absolutely necessary...
-- beautifyVariance :: Double -> Double -> Double
-- beautifyVariance _ x
--   | x < 0 = error "beautifyVariance: Variance is negative."
--   | x < eps = eps
--   | otherwise = x
--   where
--     eps = 1e-4

-- -- Analyze the covariance matrix and change problematic values. This step is
-- -- awful but necessary when there are not enough samples from the posterior.
-- beautifyCovarianceMatrix :: L.Matrix Double -> L.Matrix Double
-- beautifyCovarianceMatrix m =
--   L.accum m beautifyVariance [((i, i), 0) | i <- [0 .. nRows - 1]]
--   where
--     nRows = L.rows m

-- Read trees and extract branch lengths.
prepare :: IO ()
prepare = do
  putStrLn "Read trees."
  trsAll <- someTrees Standard fnInTrees
  let nTrees = length trsAll
  putStrLn $ show nTrees ++ " read."

  let nBurnInTrees = nTrees `div` 10
      trs = drop nBurnInTrees trsAll
  putStrLn $ "Skip a burn in of " ++ show nBurnInTrees ++ " trees."

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
          setBranches (map (either error id . toLength) $ V.toList means) (head trs)
      lvs = leaves meanTree
      trOutgroup = either error id $ outgroup (S.singleton $ head lvs) "root" meanTree
      tr = either error id $ midpoint trOutgroup
  putStrLn "The tree with mean branch lengths rooted at the midpoint:"
  print $ toNewick $ measurableToPhyloTree tr
  putStrLn $ "Save the mean tree to " <> fnMeanTree <> "."
  BL.writeFile fnMeanTree (toNewick $ measurableToPhyloTree tr)

  -- Rooting is only necessary if the obtained trees are unrooted. Here, the
  -- midpoint of the mean tree is used as root. If available, an outgroup could
  -- also be used.
  putStrLn "Root the trees at the midpoint of the mean tree."
  let outgrp = fst $ fromBipartition $ either error id $ bipartition tr
      trsRooted = map (either error id . outgroup outgrp "root") trs

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

  -- putStrLn "Beautify covariance matrix. Ouch!"
  -- let sigma = beautifyCovarianceMatrix sigmaBare
  --     variances = L.takeDiag sigma
  -- putStrLn $ "Minimum variance: " ++ show (L.minElement variances)
  -- putStrLn $ "Maximum variance: " ++ show (L.maxElement variances)
  let sigma = sigmaBare

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

  -- Use the mean tree, and the posterior means and covariances to initialize
  -- various objects.
  let -- Calibrations.
      cb = getCalibrations meanTree
      -- Constraints.
      cs = getConstraints meanTree
      -- Prior function.
      pr' = priorDistribution cb cs
      -- Likelihood function.
      lh' = likelihoodFunction mu sigmaInv logSigmaDet
      -- Proposal cycle.
      ccl' = proposals meanTree
      -- Monitor.
      mon' = monitor cb cs
      -- Starting state.
      start' = initWith meanTree

  -- Create a seed value for the random number generator. Actually, the
  -- 'create' function is deterministic, but useful during development. For
  -- real analyses, use 'createSystemRandom'.
  g <- create

  -- Construct the Markov chain.
  let s = Settings bnAnalysis burnInSpec nIterations Overwrite (SaveWithTrace 1) Debug
      c =
        -- Have a look at the 'status' function to understand the
        -- different parameters.
        chain
          pr'
          lh'
          ccl'
          mon'
          start'
          g

  -- Run the Markov chain.
  void $ mcmcWith s (MHG c)

continueMetropolisHastings :: Int -> IO ()
continueMetropolisHastings n = do
  -- Read the mean tree and the posterior means and covariances.
  meanTree <- getMeanTree
  (mu, sigmaInv, logSigmaDet) <- getData
  putStrLn "The posterior means of the branch lengths are:"
  print mu

  -- Use the mean tree, and the posterior means and covariances to initialize
  -- various objects.
  let -- Calibrations.
      cb = getCalibrations meanTree
      -- Constraints.
      cs = getConstraints meanTree
      -- Prior function.
      pr' = priorDistribution cb cs
      -- Likelihood function.
      lh' = likelihoodFunction mu sigmaInv logSigmaDet
      -- Proposal cycle.
      ccl' = proposals meanTree
      -- Monitor.
      mon' = monitor cb cs

  -- Load the MCMC status.
  s <- loadSettings (bnAnalysis ++ ".settings")
  let i = iterations s
      s' = s {iterations = i + n, executionMode = Continue}
  c <-
    loadChainWith
      pr'
      lh'
      ccl'
      mon'
      (bnAnalysis ++ ".chain")
  void $ mcmcWith s' (MHG c)

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
    -- Print usage instructions if none of the previous commands was entered.
    _ -> putStrLn "Use one command of: [prepare|run|continue N]."
