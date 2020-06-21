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
import qualified Data.Text.Lazy as T
import Lens.Micro
import Mcmc
import Numeric.Log
import Statistics.Distribution hiding (Mean)
import Statistics.Distribution.Normal
import System.Random.MWC

-- Use Int node labels for now.
type Node = Int

-- The sate space is a tree with branch lengths measured in 'Double', and 'Int'
-- node labels. We use a phantom type "a" to denote if the tree stores the
-- actual length, the mean, or the standard deviation of the branch.
-- type Tree a = AdjacencyMap (Distance Double) Node
type Tree a = AdjacencyMap Double Node

instance ToJSON (AdjacencyMap Double Node)

-- Different types of branch lengths.
type Length = Double

type Mean = Double

type StdDev = Double

-- Tree with posterior means as branches.
type MTree = Tree Mean

-- Tree with posterior standard deviations as branches.
type STree = Tree StdDev

-- The tree with current branch lengths, which is also the state of the chain,
-- this is also the used state space.
type LTree = Tree Length

getLens :: Node -> Node -> Lens' LTree Double
getLens x y = lens (g x y) (s x y)
  where
    g = edgeLabel
    s n m gr e = replaceEdge e n m gr

-- For now, use a completely uninformative prior.
prior :: a -> Log Double
prior = const 1

llhBranch :: Mean -> StdDev -> Length -> Log Double
llhBranch m v l
  | l <= 0 = negInf
  | otherwise = Exp $ logDensity (normalDistr m v) l

llh :: MTree -> STree -> LTree -> Log Double
llh mt vt lt = product $ zipWith3 llhBranch ms vs ls
  where
    ms = map (^. _1) $ edgeList mt
    vs = map (^. _1) $ edgeList vt
    ls = map (^. _1) $ edgeList lt

allEdges :: LTree -> [(Node, Node)]
allEdges = map (\(_, x, y) -> (x, y)) . edgeList

slideBr :: Node -> Node -> Move LTree
slideBr x y = slideStandard n 1 (getLens x y) True
  where
    n = "Slide edge " <> show (x, y)

moveCycle :: Cycle LTree
moveCycle = fromList $ map (uncurry slideBr) (allEdges lTree)

-- A stupid Dioid instance for 'Double'. Especially 'mempty' is problematic.
instance Semigroup Double where
  (<>) = min

instance Monoid Double where
  mempty = 1 / 0

instance Semiring Double where
  one = 0

  (<.>) = (+)

instance Dioid Double

lTree :: LTree
lTree = (0 -< 1.0 >- 1) + (0 -< 2.0 >- 2)

mTree :: MTree
mTree = (0 -< 5.0 >- 1) + (0 -< 10.0 >- 2)

vTree :: MTree
vTree = (0 -< 2.0 >- 1) + (0 -< 2.0 >- 2)

mons :: [MonitorParameter LTree]
mons = [monitorRealFloat (n x y) (getLens x y) | (x, y) <- allEdges lTree]
  where
    n x y = T.pack $ show (x, y)

monStd :: MonitorStdOut LTree
monStd = monitorStdOut mons 100

monFile :: MonitorFile LTree
monFile = monitorFile "Branches" mons 10

monBs :: [MonitorParameterBatch LTree]
monBs =
  [monitorBatchMeanRealFloat (n x y) (getLens x y) | (x, y) <- allEdges lTree]
  where
    n x y = T.pack $ "Mean " <> show (x, y)

monBatch :: MonitorBatch LTree
monBatch = monitorBatch "Branches" monBs 200

mon :: Monitor LTree
mon = Monitor monStd [monFile] [monBatch]

nBurn :: Maybe Int
nBurn = Just 2000

nIter :: Int
nIter = 10000

main :: IO ()
main = do
  g <- create
  let s =
        status
          "ApproximatePhylogeneticLikelihood"
          prior
          (llh mTree vTree)
          moveCycle
          mon
          lTree
          nBurn
          Nothing
          nIter
          g
  void $ mh s
