{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE RankNTypes #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Main
-- Description :  Approximate phylogenetic likelihood using normal distribution
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Jun 10 22:07:11 2020.
module Main
  ( main,
  )
where

import Algebra.Graph.Label
import Algebra.Graph.Labelled.AdjacencyMap
import Control.Monad
import Data.Aeson
import Data.Maybe
import Lens.Micro
import Lens.Micro.Extras
import Mcmc
import Statistics.Distribution hiding (Mean)
import Statistics.Distribution.Normal
import System.Environment
import System.Random.Stateful

-- Node (vertex) label.
type Node = Int

-- Branch length.
type Length = Double

-- The state space is a tree. We apply the branch length to the 'Distance'
-- data constructor. In this way, we tell the graph algorithms that
--
-- - branches lengths have to be added together when calculating the length of a
--   path on the graph;
--
-- - the shortest branch has to be taken if more branches connect the same set
--   of two nodes.
--
-- There are other edge label types, such as 'Capacity' (e.g., the capacity of a
-- road from one Node to another).
--
-- We use a phantom type @a@ to denote that the tree stores actual branch
-- lengths, or the means and standard deviations calculated previously.
type Tree a = AdjacencyMap (Distance Length) Node

-- We want to save and restore the Markov chain. This is done by converting
-- trees to and from JSON objects. It is unfortunate, that we have to tell
-- Haskell how to convert 'Tree' as well as 'Distance' objects.
instance ToJSON (Tree a)

instance FromJSON (Tree a)

instance ToJSON (Distance Double) where
  toJSON = toJSON . fromD
  toEncoding = toEncoding . fromD

instance FromJSON (Distance Double) where
  parseJSON d = toD <$> parseJSON d

-- Branch accessor functions (lenses).
getLens :: Node -> Node -> Lens' (Tree a) Double
getLens x y = lens (g x y) (s x y)
  where
    g v w = fromD . edgeLabel v w
    s n m gr e = replaceEdge (toD e) n m gr

-- Uninformative prior.
pr :: PriorFunction a
pr = const 1

-- Branch likelihood, a normal distribution with given mean and standard deviation.
lhBranch :: Mean Double -> StandardDeviation Double -> Length -> Likelihood
lhBranch m s l
  | l <= 0 = 0
  | otherwise = Exp $ logDensity (normalDistr m s) l

-- Likelihood of the tree for two given trees containing the branch length
-- means and standard deviations.
lh :: Tree (Mean Double) -> Tree (StandardDeviation Double) -> Tree Length -> Likelihood
lh mt st lt = product $ zipWith3 lhBranch (getEdgeLabels mt) (getEdgeLabels st) (getEdgeLabels lt)

-- Get a list of edge labels of a tree.
getEdgeLabels :: Tree a -> [Double]
getEdgeLabels = map (fromD . (^. _1)) . edgeList

-- Get a list of edges without labels of a tree.
getEdges :: Tree a -> [(Node, Node)]
getEdges = map (\(_, x, y) -> (x, y)) . edgeList

-- Sliding proposal changing the length of a branch.
slideBranch :: Node -> Node -> Proposal (Tree Length)
slideBranch x y = getLens x y @~ slideSymmetric 1.0 n (pWeight 1) Tune
  where
    n = PName $ "Edge " <> show (x, y)

-- Bactrian proposal changing the length of a branch.
bactrianBranch :: Node -> Node -> Proposal (Tree Length)
bactrianBranch x y = getLens x y @~ slideBactrian 0.9 1.0 n (pWeight 1) Tune
  where
    n = PName $ "Edge " <> show (x, y)

-- Collect all sliding proposals into a cycle.
cc :: Cycle (Tree Length)
cc =
  cycleFromList $
    map (uncurry slideBranch) (getEdges startingTree)
      ++ map (uncurry bactrianBranch) (getEdges startingTree)

-- Convert a number to a distance. Assume that the given number is positive
-- and finite.
toD :: (Num a) => a -> Distance a
toD = distance . unsafeFinite

-- Extract the number from a distance. Assume that the distance is finite.
fromD :: Distance a -> a
fromD = fromJust . getFinite . getDistance

-- Starting tree; let's start far away from the truth.
startingTree :: Tree Length
startingTree = (0 -< toD 1.0 >- 1) + (0 -< toD 2.0 >- 2)

-- Tree storing the (true) posterior means.
meanTree :: Tree (Mean Double)
meanTree = (0 -< toD 5.0 >- 1) + (0 -< toD 10.0 >- 2)

-- Tree storing the (true) posterior standard deviations.
stdDevTree :: Tree (StandardDeviation Double)
stdDevTree = (0 -< toD 2.0 >- 1) + (0 -< toD 2.0 >- 2)

-- Branch length monitors.
branchMons :: [MonitorParameter (Tree Length)]
branchMons =
  [ view (getLens x y) >$< monitorDouble (n x y)
    | (x, y) <- getEdges startingTree
  ]
  where
    n x y = show (x, y)

-- Monitor branch lengths to standard output.
monStdOut :: MonitorStdOut (Tree Length)
monStdOut = monitorStdOut branchMons 100

-- Monitor branch lengths to a file.
monFile :: MonitorFile (Tree Length)
monFile = monitorFile "branches" branchMons 10

-- Monitor batch means of branch lengths.
branchBatchMons :: [MonitorParameterBatch (Tree Length)]
branchBatchMons =
  [view (getLens x y) >$< monitorBatchMean (n x y) | (x, y) <- getEdges startingTree]
  where
    n x y = "Mean " <> show (x, y)

-- Monitor batch means of branch lengths to a file.
monBatch :: MonitorBatch (Tree Length)
monBatch = monitorBatch "branches" branchBatchMons 100

-- Combine the monitors.
mon :: Monitor (Tree Length)
mon = Monitor monStdOut [monFile] [monBatch]

analysisName :: AnalysisName
analysisName = AnalysisName "plh"

-- The main program.
main :: IO ()
main = do
  -- Get the command line arguments and check if a new chain is to be started,
  -- or a previously finished, and saved chain is to be continued.
  args <- getArgs
  case args of
    [] -> do
      let g = mkStdGen 0
      -- Combine all the objects defined above.
      let s =
            Settings
              analysisName
              (BurnInWithAutoTuning 4000 200)
              (Iterations 20000)
              TraceAuto
              Overwrite
              Sequential
              Save
              LogStdOutAndFile
              Info
      -- Initialize the Metropolis-Hastings-Green algorithm.
      a <- mhg s pr (lh meanTree stdDevTree) cc mon startingTree g
      -- Run the MCMC sampler.
      void $ mcmc s a
    ["continue", nStr] -> do
      -- Load a previously finished, and saved chain. We have to give the prior
      -- and likelihood functions, as well as the proposals and the monitors,
      -- because those cannot be saved consistently. The loading function checks
      -- the values of the prior and likelihood functions of the last state, to
      -- minimize the chance that wrong functions are used by accident. It also
      -- sets the tuning parameters of the proposals in the cycle. Using different
      -- proposals in the cycle, or using different monitors may lead to undefined
      -- behavior and is not supported.
      s <- settingsLoad analysisName
      a <- mhgLoad pr (lh meanTree stdDevTree) cc mon analysisName
      -- Continue the MCMC sampler for the given number of iterations.
      void $ mcmcContinue (read nStr) s a
    xs -> do
      p <- getProgName
      putStrLn $ "usage: " ++ p
      putStrLn $ "   or: " ++ p ++ " continue N"
      error $ "could not read command line arguments: " ++ unwords xs
