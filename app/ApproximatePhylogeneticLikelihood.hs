{- |
Module      :  ApproximatePhylogeneticLikelihood
Description :  Approximate phylogenetic likelihood
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed Jun 10 22:07:11 2020.

-}

module ApproximatePhylogeneticLikelihood
  ()
where

import           Algebra.Graph.Label
import           Algebra.Graph.Labelled
import           Data.Maybe
import           Lens.Micro
import           Numeric.Log
import           Statistics.Distribution hiding ( Mean )
import           Statistics.Distribution.Normal

import           Mcmc.Move
import           Mcmc.Move.Generic

type Length = Double
type Mean = Double
type StdDev = Double

-- Use Int node labels for now.
type Node = Int

-- Tree with posterior means as branches.
type MTree = Graph (Distance Mean) Node

-- Tree with posterior standard deviations as branches.
type VTree = Graph (Distance StdDev) Node

-- The tree with current branch lengths, which is also the state of the chain.
type LTree = Graph (Distance Length) Node
type I = LTree

-- For now, use a completely uninformative prior.
prior :: a -> Log Double
prior = const 1

negInf :: Fractional a => a
negInf = -(1 / 0)

llhBranch :: Mean -> StdDev -> Length -> Log Double
llhBranch m v l | l <= 0    = Exp negInf
                | otherwise = Exp $ logDensity (normalDistr m v) l

llh :: MTree -> VTree -> I -> Log Double
llh mt vt lt = product $ zipWith3 llhBranch ms vs ls
 where
  ms = map (fromD . (^. _1)) $ edgeList mt
  vs = map (fromD . (^. _1)) $ edgeList vt
  ls = map (fromD . (^. _1)) $ edgeList lt

setter :: (Monoid e, Eq e, Ord a) => a -> a -> e -> Graph e a -> Graph e a
setter x y e = replaceEdge e x y

getter :: (Monoid e, Eq a) => a -> a -> Graph e a -> e
getter = edgeLabel

allEdges :: (Eq e, Monoid e, Ord a) => Graph e a -> [(a, a)]
allEdges = map (\(_, x, y) -> (x, y)) . edgeList

slide :: Node -> Node -> Move I
slide x y = Move n 1 slideSimple Nothing
 where
  n           = show (x, y)
  slideSimple = moveGenericContinuous (setter x y . distance . unsafeFinite)
                                      (fromD . getter x y)
                                      (normalDistr 0 1.0)
                                      (+)
                                      (-)

moveCycle :: Cycle I
moveCycle = fromList $ map (uncurry slide) (allEdges lTree)

toD :: a -> Distance a
toD = distance . unsafeFinite

fromD :: Distance a -> a
fromD = fromJust. getFinite . getDistance

lTree :: LTree
lTree = (0 -< toD 1.0 >- 1) + (0 -< toD 2.0 >- 2)

mTree :: MTree
mTree = (0 -< toD 5.0 >- 1) + (0 -< toD 10.0 >- 2)

vTree :: MTree
vTree = (0 -< toD 2.0 >- 1) + (0 -< toD 2.0 >- 2)

main :: IO ()
main = do
  g           <- create
  mu_observed <- arrowMean g
  -- putStrLn $ "True parameter: " <> show mu_observed
  let s = status 0.01 prior (likelihood mu_observed) moveCycle mon g
  void $ mh nBurn nAutoTune nIter s
