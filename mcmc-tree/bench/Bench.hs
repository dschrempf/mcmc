{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Bench
-- Description :  Benchmark common functions
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Nov  3 14:47:08 2020.
module Main
  ( main,
  )
where

import Control.Lens
import Criterion.Main
import Data.Maybe
import ELynx.Tree
import Mcmc.Tree

--   let bf1 =
--         toTree . insertLabel "Bla"
--           . fromMaybe (error "Path does not lead to a leaf.")
--           . goPath pth
--           . fromTree
--   putStrLn $ "Change a leaf: " <> show (bf1 tr) <> "."
--   putStrLn "Benchmark change a leaf."
--   benchmark $ nf bf1 tr
--   let bf2 =
--         label . current
--           . fromMaybe (error "Path does not lead to a leaf.")
--           . goPath pth
--           . fromTree
--   putStrLn $ "Leaf to get: " <> show (bf2 tr) <> "."
--   putStrLn "Benchmark get a leaf."
--   benchmark $ nf bf2 tr
--   putStrLn "Benchmark calculation of prior."
--   let i = initWith tr
--       pr' = priorFunction (getCalibrations tr) (getConstraints tr)
--   benchmark $ nf pr' i
--   (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' fnData
--   let sigmaInv = L.fromRows sigmaInvRows
--       lh' = likelihoodFunction mu sigmaInv logSigmaDet
--   putStrLn "Benchmark calculation of likelihood."
--   benchmark $ nf lh' i

--   putStrLn "Benchmark identify."
--   benchmark $ nf identify tr

changeLeaf :: Lens' (Tree e a) a -> a -> Tree e a -> Tree e a
changeLeaf l x t = t & l .~ x

getPath :: (Show a, Eq a) => a -> Tree e a -> Path
getPath x t = fst $ fromMaybe err $ ifind (\_ n -> n == x) t
  where
    err = error $ show x <> " not found."

main :: IO ()
main = do
  tr <- oneTree Standard "data/bench.tree"
  let nameGn = "Gnetum_montanum"
      pthGn = getPath nameGn tr
  let nameBr = "Brachypodium_distachyon"
      pthBr = getPath nameBr tr
  -- -- Some debugging.
  -- putStrLn $ "The path to \"Gnetum_montanum\" is: " <> show pth <> "."
  -- print $ toNewick $ measurableToPhyloTree tr
  -- let tr' = changeLeaf (subTreeAtL pth . root) "" tr
  -- print $ toNewick $ measurableToPhyloTree tr'
  defaultMain
    -- Optimized lenses are around 10 percent faster for this tree.
    [ bgroup
        "lens"
        [ bench "change leaf Gn, optimized lens" $
            nf (changeLeaf (subTreeAtL pthGn . labelL) "") tr,
          bench "change leaf Br, optimized lens" $
            nf (changeLeaf (subTreeAtL pthBr . labelL) "") tr
        ]
    ]
-- Post scriptum:
--
-- Benchmarks with criterion indicate the following:
--
-- 1. Algebraic graphs are slow, but the adjacency map with 'Int' node labels
-- and branch labels is reasonably fast.
--
-- 2. Lenses are fast, but don't allow accessing children, or the parent.
--
-- 3. Zippers are very slow.
--
-- * Lenses:
--
-- benchmarking lens/change leaf Gn, optimized lens
-- time                 5.870 μs   (5.852 μs .. 5.906 μs)
--                      1.000 R²   (1.000 R² .. 1.000 R²)
-- mean                 5.864 μs   (5.857 μs .. 5.889 μs)
-- std dev              42.95 ns   (15.43 ns .. 86.41 ns)
--
-- benchmarking lens/change leaf Br, optimized lens
-- time                 6.285 μs   (6.228 μs .. 6.365 μs)
--                      0.999 R²   (0.997 R² .. 1.000 R²)
-- mean                 6.309 μs   (6.227 μs .. 6.534 μs)
-- std dev              415.1 ns   (155.0 ns .. 820.7 ns)
-- variance introduced by outliers: 74% (severely inflated)
--
-- * Adjacency maps
--
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
-- * Zippers
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
