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

-- import AdjacencyIntMap
-- import Data.List
-- import Lens.Micro.Platform
-- import Mcmc
-- import Numeric.LinearAlgebra (Matrix, (<#), (<.>))
-- import System.Random.MWC

import Control.Lens
import Control.Monad
import Criterion
import Data.Aeson
import Data.Bifunctor
import Data.ByteString.Lazy.Char8 (ByteString)
import qualified Data.ByteString.Lazy.Char8 as L
import Data.List
import Data.Maybe
import Data.Set (Set)
import qualified Data.Set as S
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import ELynx.Data.Tree
import qualified ELynx.Data.Topology.Rooted as T
import ELynx.Export.Tree.Newick
import GHC.Generics
import Numeric.LinearAlgebra (Matrix)
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
    timeTree :: Tree Double Int,
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

-- Prior.
pr :: I -> Log Double
pr (I l h t k m r) =
  product'
    [ exponentialWith 1.0 l,
      exponentialWith 10.0 h,
      branchesWith (exponentialWith l) t,
      exponentialWith 10.0 k,
      gammaWith k (1/k) m,
      branchesWith (gammaWith k (1/k)) r
    ]

instance ToJSON I
instance FromJSON I


fn :: FilePath
fn = "mcmc-examples/PhylogeneticLikelihoodMultivariate/data/plants_1.treelist.gz"

outgroups :: Set ByteString
outgroups =
  S.fromList
    [ "Nu_advena",
      "Gi_biloba",
      "Cy_micholi",
      "Ps_nudum",
      "Sc_dissect",
      "Op_vulgatu",
      "Pi_radiata",
      "Pi_pondero",
      "Pi_jeffrey",
      "Ce_libani",
      "Th_elegans",
      "Cy_spinulo",
      "Eq_diffusu",
      "Ar_thalian",
      "Yu_filamen",
      "Po_acrosti",
      "Di_truncat",
      "Th_acumina",
      "Ho_pycnoca",
      "Gy_dryopte",
      "Cy_utahens",
      "Cy_fragili",
      "Pt_vittata",
      "Pt_ensigor",
      "Ad_tenerum",
      "Ad_aleutic",
      "Di_villosa",
      "Pe_borboni",
      "Po_trichoc",
      "Po_euphrat",
      "Ep_sinica",
      "Sc_vertici",
      "Pi_parvifl",
      "Ze_mays",
      "So_bicolor",
      "Ac_america",
      "Sm_bona",
      "Co_autumna",
      "Lu_polyphy",
      "Oe_specios",
      "Oe_rosea",
      "La_trident",
      "Ne_nucifer",
      "Vi_vinifer",
      "Ko_scopari",
      "Lu_angusti",
      "Ro_chinens",
      "Bo_nivea",
      "Hi_cannabi",
      "Ca_papaya",
      "Di_malabar",
      "Ta_parthen",
      "In_heleniu",
      "Ro_officin",
      "So_tuberos",
      "Ip_purpure",
      "Ca_roseus",
      "Al_cathart",
      "Ju_scopulo",
      "Cu_lanceol",
      "Ho_cordata",
      "Sa_bermuda",
      "Or_sativa",
      "Br_distach",
      "Zo_marina",
      "Po_peltatu",
      "Es_califor",
      "Ka_heteroc",
      "Il_parvifl",
      "Ma_attenua",
      "Da_nodosa",
      "Di_conjuga",
      "Ta_baccata",
      "Sa_henryi",
      "Am_trichop",
      "Pr_andina",
      "Gn_montanu",
      "Il_florida",
      "Eq_hymale",
      "Sa_glabra"
    ]

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

getPosteriorMatrix :: [Tree Double a] -> Matrix Double
getPosteriorMatrix = L.fromRows . map (sumFirstTwo . getBranches)

-- -- Phylogenetic likelihood using a multivariate normal distribution. See
-- -- https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
-- --
-- -- The constant @k * log (2*pi)@ was left out on purpose.
-- --
-- -- lh meanVector invertedCovarianceMatrix logOfDeterminantOfCovarianceMatrix
-- lh :: Vector Double -> Matrix Double -> Double -> I -> Log Double
-- lh mu sigmaInv logSigmaDet xs = Exp $ (-0.5) * (logSigmaDet + ((dxs <# sigmaInv) <.> dxs))
--   where
--     dxs = xs - mu

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
    ["inspect"] ->
      do
        tr <- oneTree fn
        putStrLn $ "The tree has " <> show (length $ leaves tr) <> " leaves."
        let trRooted = either error id $ outgroup outgroups "root" tr
        L.putStrLn $ toNewick $ lengthToPhyloTree trRooted
        print $ getBranches $ first fromLength trRooted
        let xs = itoList trRooted
        print xs
        let (pth, _) = fromMaybe (error "Gn_montanu not found.") $ ifind (\_ n -> n == "Gn_montanu") trRooted
        print pth
        let bf1 =
              toTree . insertLabel "BLAAAAAAAAAAAAAAA"
                . fromMaybe (error "Dohh")
                . goPath pth
                . fromTree
        let bf2 =
              label . current
                . fromMaybe (error "Dohh")
                . goPath pth
                . fromTree
        print $ bf1 trRooted
        print $ bf2 trRooted
        benchmark $ nf bf1 trRooted
        benchmark $ nf bf2 trRooted
    ["read"] -> do
      putStrLn "Read trees."
      trs <- someTrees fn
      let l = length $ nub $ map T.fromLabeledTree trs
      unless (l == 1) (error "Trees have different topologies.")
      let trsRooted = map (either error id . outgroup outgroups "root") trs

      putStrLn "Get the posterior means and the posterior covariance matrix."
      let pm = getPosteriorMatrix $ map (first fromLength) trsRooted
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
