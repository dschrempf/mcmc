{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE RankNTypes #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Main
-- Description :  Approximate phylogenetic likelihood
-- Copyright   :  (c) Dominik Schrempf, 2020
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
import Data.Aeson.Types
import Data.Maybe
import Lens.Micro
import Mcmc
import Numeric.Log
import Statistics.Distribution hiding (Mean)
import Statistics.Distribution.Normal
import System.Environment
import System.Random.MWC

-- Node (vertex) label.
type Node = Int

-- Branch length.
type Length = Double

-- Mean.
type Mean = Double

-- Standard deviation.
type StdDev = Double

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
  parseJSON d = toD <$> (parseJSON d :: Parser Double)

-- Branch accessor functions (lenses).
getLens :: Node -> Node -> Lens' (Tree a) Double
getLens x y = lens (g x y) (s x y)
  where
    g v w = fromD . edgeLabel v w
    s n m gr e = replaceEdge (toD e) n m gr

-- Uninformative prior.
pr :: a -> Log Double
pr = const 1

-- Branch likelihood, a normal distribution with given mean and standard deviation.
lhBranch :: Mean -> StdDev -> Length -> Log Double
lhBranch m s l
  | l <= 0 = pzero
  | otherwise = Exp $ logDensity (normalDistr m s) l

-- Likelihood of the tree for two given trees containing the branch length
-- means and standard deviations.
lh :: Tree Mean -> Tree StdDev -> Tree Length -> Log Double
lh mt st lt = product $ zipWith3 lhBranch (getEdgeLabels mt) (getEdgeLabels st) (getEdgeLabels lt)

-- Get a list of edge labels of a tree.
getEdgeLabels :: Tree a -> [Double]
getEdgeLabels = map (fromD . (^. _1)) . edgeList

-- Get a list of edges without labels of a tree.
getEdges :: Tree a -> [(Node, Node)]
getEdges = map (\(_, x, y) -> (x, y)) . edgeList

-- Sliding move changing the length of a branch.
slideBranch :: Node -> Node -> Move (Tree Length)
slideBranch x y = slideSymmetric n 1 (getLens x y) 1.0 True
  where
    n = "Slide edge " <> show (x, y)

-- Bactrian move changing the length of a branch.
bactrianBranch :: Node -> Node -> Move (Tree Length)
bactrianBranch x y = slideBactrian n 1 (getLens x y) 0.9 1.0 True
  where
    n = "Bactrian edge " <> show (x, y)

-- Collect all sliding moves into a cycle.
moveCycle :: Cycle (Tree Length)
moveCycle =
  fromList $
    map (uncurry slideBranch) (getEdges startingTree)
      ++ map (uncurry bactrianBranch) (getEdges startingTree)

-- Convert a number to a distance. Assume that the given number is positive
-- and finite.
toD :: Num a => a -> Distance a
toD = distance . unsafeFinite

-- Extract the number from a distance. Assume that the distance is finite.
fromD :: Distance a -> a
fromD = fromJust . getFinite . getDistance

-- Starting tree; let's start far away from the truth.
startingTree :: Tree Length
startingTree = (0 -< toD 1.0 >- 1) + (0 -< toD 2.0 >- 2)

-- Tree storing the (true) posterior means.
meanTree :: Tree Mean
meanTree = (0 -< toD 5.0 >- 1) + (0 -< toD 10.0 >- 2)

-- Tree storing the (true) posterior standard deviations.
stdDevTree :: Tree StdDev
stdDevTree = (0 -< toD 2.0 >- 1) + (0 -< toD 2.0 >- 2)

-- Branch length monitors.
branchMons :: [MonitorParameter (Tree Length)]
branchMons = [monitorRealFloat (n x y) (getLens x y) | (x, y) <- getEdges startingTree]
  where
    n x y = show (x, y)

-- Monitor branch lengths to standard output.
monStdOut :: MonitorStdOut (Tree Length)
monStdOut = monitorStdOut branchMons 100

-- Monitor branch lengths to a file.
monFile :: MonitorFile (Tree Length)
monFile = monitorFile "Branches" branchMons 10

-- Monitor batch means of branch lengths.
branchBatchMons :: [MonitorParameterBatch (Tree Length)]
branchBatchMons =
  [monitorBatchMeanRealFloat (n x y) (getLens x y) | (x, y) <- getEdges startingTree]
  where
    n x y = "Mean " <> show (x, y)

-- Monitor batch means of branch lengths to a file.
monBatch :: MonitorBatch (Tree Length)
monBatch = monitorBatch "Branches" branchBatchMons 100

-- Combine the monitors.
mon :: Monitor (Tree Length)
mon = Monitor monStdOut [monFile] [monBatch]

-- Number of burn in iterations.
nBurnIn :: Maybe Int
nBurnIn = Just 4000

-- Auto tuning period.
nAutoTune :: Maybe Int
nAutoTune = Just 200

-- Number of Metropolis-Hasting iterations after burn in.
nIterations :: Int
nIterations = 20000

-- Name of the chain; used as prefix of the output file names.
nm :: String
nm = "ApproximatePhylogeneticLikelihood"

-- The main program.
main :: IO ()
main = do
  -- Get the command line arguments and check if a new chain is to be started,
  -- or a previously finished, and saved chain is to be continued.
  args <- getArgs
  case args of
    [] -> do
      g <- create
      -- Combine all the objects defined above.
      let s = status nm pr (lh meanTree stdDevTree) moveCycle mon startingTree nBurnIn nAutoTune nIterations g
      -- Run the Markov chain Monte Carlo sampler using the Metropolis-Hastings algorithm.
      void $ mh s
    ["continue", nStr] -> do
      -- Load a previously finished, and saved chain. We have to give the prior
      -- and likelihood functions, as well as the move cycle and the monitors,
      -- because those cannot be saved consistently. The loading function checks
      -- the values of the prior and likelihood functions of the last state, to
      -- minimize the chance that wrong functions are used by accident. It also
      -- sets the tuning parameters of the moves in the cycle. Using different
      -- moves in the cycle, or using different monitors may lead to undefined
      -- behavior and is not supported.
      s <- loadStatus pr (lh meanTree stdDevTree) moveCycle mon (nm ++ ".mcmc")
      -- Continue the chain for the given number of iterations.
      void $ mhContinue (read nStr) s
    xs -> do
      p <- getProgName
      putStrLn $ "usage: " ++ p
      putStrLn $ "   or: " ++ p ++ " continue N"
      error $ "could not read command line arguments: " ++ unwords xs
