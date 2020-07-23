{-# LANGUAGE DeriveGeneric #-}
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
-- The trees are read from a data file which is given relative to the @mcmc@ git
-- repository base directory. Hence, the compiled binary has to be executed from
-- this directory.
module Main
  ( main,
  )
where

import Control.Lens hiding ((<.>))
import Control.Monad
import Criterion
import Data.Aeson
import Data.Bifunctor
import Data.Bitraversable
import qualified Data.ByteString.Lazy.Char8 as L
import Data.List
import Data.Maybe
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import Debug.Trace
import qualified ELynx.Data.Topology.Rooted as T
import ELynx.Data.Tree
import ELynx.Export.Tree.Newick
import GHC.Generics
import Numeric.LinearAlgebra (Matrix, (<#), (<.>))
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import Prior
import System.Environment
import Tree

-- State space.
data I = I
  { -- Birth rate parameter of time tree.
    timeBirthRate :: Double,
    -- Height of root node measured in units of time.
    timeRootHeight :: Double,
    -- Time tree.
    -- TODO: The types are all wrong.
    timeTree :: Tree Length Int,
    -- Shape parameter k of gamma distribution of rate parameters. The scale
    -- parameter is determined such that the mean of the gamma distribution is
    -- 1.
    rateGammaShape :: Double,
    -- Global mean of rate parameters.
    rateMean :: Double,
    -- Rate tree.
    rateTree :: Tree Double Int
  }
  deriving (Generic)

instance ToJSON I

instance FromJSON I

-- Initial state.
initWith :: Tree Length Int -> I
initWith t =
  I
    { timeBirthRate = 1.0,
      timeRootHeight = height t',
      timeTree = t',
      rateGammaShape = 1.0,
      rateMean = 1.0,
      rateTree = first (const 1.0) t
    }
  where
    t' = makeUltrametric t

-- Prior.
pr :: I -> Log Double
pr (I l h t k m r) =
  product'
    [ exponentialWith 1.0 l,
      exponentialWith 10.0 h,
      branchesWith (exponentialWith l) (first fromLength t),
      exponentialWith 10.0 k,
      gammaWith k (1 / k) m,
      branchesWith (gammaWith k (1 / k)) r
    ]

fnTreeList :: FilePath
fnTreeList = "mcmc-examples/PhylogeneticLikelihoodMultivariate/data/plants_1.treelist.gz"

fnMeanTree :: FilePath
fnMeanTree = "mcmc-examples/PhylogeneticLikelihoodMultivariate/data/Mean.tree"

-- Get all branches of the tree such that the two branches leading to the root
-- are the first two entries of the vector. Ignore the root branch.
getBranches :: Tree Double a -> Vector Double
getBranches (Node _ _ [l, r]) = V.fromList $ head ls : head rs : tail ls ++ tail rs
  where
    ls = branches l
    rs = branches r
getBranches _ = error "getBranches: Root node is not bifurcating."

-- Sum the first two elements of a vector.
sumFirstTwo :: Vector Double -> Vector Double
sumFirstTwo v = (v V.! 0 + v V.! 1) `V.cons` V.drop 2 v

getPosteriorMatrixRooted :: [Tree Double a] -> Matrix Double
getPosteriorMatrixRooted = L.fromRows . map (sumFirstTwo . getBranches)

getPosteriorMatrix :: [Tree Double a] -> Matrix Double
getPosteriorMatrix = L.fromRows . map (V.fromList . branches)

-- | Apply a function with different effect on each node to a 'Traversable'.
-- Based on https://stackoverflow.com/a/41523456.
tZipWith :: Bitraversable t => [a] -> t b c -> Maybe (t a c)
tZipWith xs = bisequenceA . snd . bimapAccumL pair noChange xs
  where
    pair [] _ = ([], Nothing)
    pair (y : ys) _ = (ys, Just y)
    noChange ys z = (ys, Just z)

-- Phylogenetic likelihood using a multivariate normal distribution. See
-- https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
--
-- The constant @k * log (2*pi)@ was left out on purpose.
--
-- lh meanVector invertedCovarianceMatrix logOfDeterminantOfCovarianceMatrix
lh :: Vector Double -> Matrix Double -> Double -> I -> Log Double
lh mu sigmaInv logSigmaDet x = traceShow dxs $ Exp $ (-0.5) * (logSigmaDet + ((dxs <# sigmaInv) <.> dxs))
  where
    times = getBranches $ first fromLength $ timeTree x
    rates = getBranches $ rateTree x
    multiplier = timeRootHeight x * rateMean x
    distances = sumFirstTwo $ V.map (* multiplier) $ V.zipWith (*) times rates
    dxs = distances - mu

-- -- Slide branch with given index.
-- slideBranch :: Int -> Proposal I
-- slideBranch i = slideSymmetric n 1 (singular $ ix i) 0.01 True
--   where
--     n = "Slide branch " <> show i

-- proposals :: I -> Cycle I
-- proposals v = fromList [slideBranch i | i <- [0 .. k]]
--   where
--     k = V.length v - 1

-- -- Branch length monitors.
-- branchMons :: I -> [MonitorParameter I]
-- branchMons v = [monitorRealFloat (n i) (singular $ ix i) | i <- [0 .. k]]
--   where
--     n i = "Branch " <> show i
--     k = V.length v - 1

-- mon :: I -> Monitor I
-- mon v = Monitor (monitorStdOut (take 3 bs) 50) [monitorFile "Branches" bs 10] []
--   where
--     bs = branchMons v

-- -- Number of burn in iterations.
-- nBurnIn :: Maybe Int
-- nBurnIn = Just 1600

-- -- Auto tuning period.
-- nAutoTune :: Maybe Int
-- nAutoTune = Just 200

-- -- Number of Metropolis-Hasting iterations after burn in.
-- nIterations :: Int
-- nIterations = 10000

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

fnData :: String
fnData = "plh-multivariate.data"

main :: IO ()
main = do
  as <- getArgs

  case as of
    ["mean"] -> do
      putStrLn "Read trees; skip a burn in of 1000 trees."
      trs <- drop 1000 <$> someTrees fnTreeList
      let l = length $ nub $ map T.fromLabeledTree trs
      unless (l == 1) (error "Trees have different topologies.")
      let pm = getPosteriorMatrix $ map (first fromLength) trs
          (mu, _) = L.meanCov pm
      print mu
      let meanTree =
            fromMaybe (error "Could not label tree with mean branch lengths") $
              tZipWith (V.toList mu) (head trs)
      putStrLn "The tree with mean branch lengths:"
      print meanTree
      putStrLn $ "Save the tree to " <> fnMeanTree <> "."
      L.writeFile fnMeanTree (toNewick $ lengthToPhyloTree $ first Length meanTree)
    ["inspect"] ->
      do
        -- TODO: Midpoint root.
        tr <- oneTree fnTreeList
        putStrLn $ "The tree has " <> show (length $ leaves tr) <> " leaves."
        -- let trRooted = either error id $ outgroup outgroups "root" tr
        putStrLn "The rooted tree is:"
        L.putStrLn $ toNewick $ lengthToPhyloTree tr
        putStrLn "The branch lengths are:"
        print $ getBranches $ first fromLength tr
        let (pth, _) =
              fromMaybe (error "Gn_montanu not found.") $
                ifind (\_ n -> n == "Gn_montanu") tr
        putStrLn "The path to \"Gn_montanu\" is:"
        print pth
        -- let bf1 =
        --       toTree . insertLabel "Bla"
        --         . fromMaybe (error "Path does not lead to a leaf.")
        --         . goPath pth
        --         . fromTree
        -- putStrLn $ "Change a leaf: " <> show (bf1 trRooted) <> "."
        -- putStrLn "Benchmark change a leaf."
        -- benchmark $ nf bf1 trRooted
        -- let bf2 =
        --       label . current
        --         . fromMaybe (error "Path does not lead to a leaf.")
        --         . goPath pth
        --         . fromTree
        -- putStrLn $ "Leaf to get: " <> show (bf2 trRooted) <> "."
        -- putStrLn "Benchmark get a leaf."
        -- benchmark $ nf bf2 trRooted
        let i = initWith $ identify tr
        putStrLn $ "Test if time tree is ultrametric: " <> show (ultrametric $ timeTree i)
        putStrLn $ "Initial prior: " <> show (pr i) <> "."
        putStrLn "Benchmark calculation of prior."
        benchmark $ nf pr i
        putStrLn "Load posterior means and covariances."
        (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
        let sigmaInv = L.fromRows sigmaInvRows
            lh' = lh mu sigmaInv logSigmaDet
        -- TODO: Likelihood is zero?
        putStrLn $ "Initial likelihood: " <> show (lh' i) <> "."
    -- putStrLn "Benchmark calculation of likelihood."
    -- benchmark $ nf lh' i
    ["read"] -> do
      tr <- oneTree fnTreeList
      let outgroups = fst $ fromBipartition $ either error id $ bipartition tr
      putStrLn "Read trees; skip a burn in of 1000 trees."
      trs <- drop 1000 <$> someTrees fnTreeList
      -- let l = length $ nub $ map T.fromLabeledTree trs
      -- unless (l == 1) (error "Trees have different topologies.")
      let trsRooted = map (either error id . outgroup outgroups "root") trs

      putStrLn "Get the posterior means and the posterior covariance matrix."
      let pm = getPosteriorMatrixRooted $ map (first fromLength) trsRooted
          (mu, sigma) = L.meanCov pm
          (sigmaInv, (logSigmaDet, _)) = L.invlndet $ L.unSym sigma

      putStrLn "The posterior means of the branch lengths are:"
      print mu

      putStrLn $ "Save posterior means an covariances to " <> fnData <> "."
      encodeFile fnData (mu, L.toRows sigmaInv, logSigmaDet)
    _ -> putStrLn "inspect|read"

-- _ -> do
--   g <- create
--   (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
--   let sigmaInv = L.fromRows sigmaInvRows
--   putStrLn "The posterior means of the branch lengths are:"
--   print mu
--   putStrLn "Choose a (bad) starting state for our chain."
--   let k = V.length mu
--       start = V.replicate k (1.0 :: Double)
--   print start
--     putStrLn "Construct status of the chain."
--     let s = force $ status "plh-multivariate" pr (lh mu sigmaInv logSigmaDet) (proposals start) (mon start) start nBurnIn nAutoTune nIterations g
--     void $ mh s
