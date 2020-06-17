{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Main
Description :  Approximate phylogenetic likelihood
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed Jun 10 22:07:11 2020.

-}

module Main
  ( main
  )
where

import           Algebra.Graph.Label
import           Algebra.Graph.Labelled.AdjacencyMap
import           Control.Monad
import           Data.Maybe
import qualified Data.Text.Lazy                as T
import           Lens.Micro
import           Numeric.Log
import           Statistics.Distribution hiding ( Mean )
import           Statistics.Distribution.Normal
import           System.Random.MWC

import           Mcmc

type Length = Double
type Mean = Double
type StdDev = Double

-- Use Int node labels for now.
type Node = Int

-- Tree with posterior means as branches.
type MTree = AdjacencyMap (Distance Mean) Node

-- Tree with posterior standard deviations as branches.
type VTree = AdjacencyMap (Distance StdDev) Node

-- The tree with current branch lengths, which is also the state of the chain.
type LTree = AdjacencyMap (Distance Length) Node
type I = LTree

getLens
  :: (Num e, Ord e, Ord a) => a -> a -> Lens' (AdjacencyMap (Distance e) a) e
getLens x y = lens (g x y) (s x y)
 where
  g n m = fromD . edgeLabel n m
  s n m gr e = replaceEdge (toD e) n m gr

-- For now, use a completely uninformative prior.
prior :: a -> Log Double
prior = const 1

llhBranch :: Mean -> StdDev -> Length -> Log Double
llhBranch m v l | l <= 0    = negInf
                | otherwise = Exp $ logDensity (normalDistr m v) l

llh :: MTree -> VTree -> I -> Log Double
llh mt vt lt = product $ zipWith3 llhBranch ms vs ls
 where
  ms = map (fromD . (^. _1)) $ edgeList mt
  vs = map (fromD . (^. _1)) $ edgeList vt
  ls = map (fromD . (^. _1)) $ edgeList lt

allEdges :: (Eq e, Monoid e, Ord a) => AdjacencyMap e a -> [(a, a)]
allEdges = map (\(_, x, y) -> (x, y)) . edgeList

slideBr :: Node -> Node -> Move I
slideBr x y = slideStandard n 1 (getLens x y) True
 where
  n           = "Slide edge " <> show (x, y)

moveCycle :: Cycle I
moveCycle = fromList $ map (uncurry slideBr) (allEdges lTree)

toD :: a -> Distance a
toD = distance . unsafeFinite

fromD :: Distance a -> a
fromD = fromJust . getFinite . getDistance

lTree :: LTree
lTree = (0 -< toD 1.0 >- 1) + (0 -< toD 2.0 >- 2)

mTree :: MTree
mTree = (0 -< toD 5.0 >- 1) + (0 -< toD 10.0 >- 2)

vTree :: MTree
vTree = (0 -< toD 2.0 >- 1) + (0 -< toD 2.0 >- 2)

mons :: [MonitorParameter I]
mons = [ monitorRealFloat (n x y) (getLens x y) | (x, y) <- allEdges lTree ]
  where n x y = T.pack $ show (x, y)

monStd :: MonitorStdOut I
monStd = monitorStdOut mons 100

monFile :: MonitorFile I
monFile = monitorFile "ApproximatePhylogeneticLikelihood.log" mons 10

mon :: Monitor I
mon = Monitor monStd [monFile] []

nBurn :: Maybe Int
nBurn = Just 2000

nIter :: Int
nIter = 10000

main :: IO ()
main = do
  g <- create
  let s = status prior (llh mTree vTree) moveCycle mon lTree g
  void $ mh nBurn Nothing nIter s
