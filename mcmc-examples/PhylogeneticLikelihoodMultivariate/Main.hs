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
-- import Numeric.Log
-- import Numeric.LinearAlgebra (Matrix, (<#), (<.>))
-- import System.Random.MWC

import Control.Lens
import Control.Monad
import Control.Zipper
import Data.Aeson
import Data.Bifoldable
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
import Numeric.LinearAlgebra (Matrix)
import qualified Numeric.LinearAlgebra as L
import System.Environment
import Tree

-- TODO: Birth-death prior on time-tree (substitution-tree)?

-- -- We condense the branch lengths into a vector.
-- data I = I
--   { -- Height of root node measured in units of time.
--     timeRoot :: Double,
--     -- Global normalization of rate parameters. XXX: Is this necessary given we
--     -- have a fully specified the gamma distribution with mean k*theta (see
--     -- below)?
--     rateMean :: Double,
--     -- First parameter k of gamma distribution of rate parameters.
--     rateGammaScale :: Double,
--     -- Second parameter theta of gamma distribution of rate parameters.
--     rateGammaShape :: Double,
--     -- Tree.
--     tree :: Tree Double Int
--   }

-- instance ToJSON I
-- instance FromJSON I

-- -- Does not work because the likelihood function will not be recomputed since
-- -- it believes the vector has not changed.
-- unsafeSet :: Int -> Vector R -> Double -> Vector R
-- unsafeSet i v x = runST $ do
--   mv <- V.unsafeThaw v
--   M.write mv i x
--   V.unsafeFreeze mv

-- -- Custom mutable lens.
-- unsafeIx :: Int -> Lens' (Vector R) R
-- unsafeIx i = lens (V.! i) (unsafeSet i)

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

-- Get the vecotr of branch lengths.
--
-- Ignore the stem, and sum the two branches leading to the root.
getBranches :: Measurable e => Tree e a -> Vector Double
getBranches (Node _ _ [l, r]) = V.fromList $ (getStem l + getStem r) : concatMap go (forest l) ++ concatMap go (forest r)
  where
    go = bifoldr' (\x acc -> getLen x : acc) (flip const) []
getBranches _ = error "getBranches: Root node is not bifurcating."

getPosteriorMatrix :: [Tree Length a] -> Matrix Double
getPosteriorMatrix = L.fromRows . map getBranches

-- -- Uniform prior. Ensuring positive branch lengths. If this is too slow, the
-- -- positiveness of branches has to be ensured by the proposals.
-- pr :: I -> Log Double
-- pr xs
--   | V.any (<= 0) xs = pzero
--   | otherwise = Exp 0

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
-- Zipper; preview node 150
-- benchmarking...
-- time                 2.513 μs   (2.351 μs .. 2.678 μs)
--                      0.978 R²   (0.972 R² .. 0.993 R²)
-- mean                 2.432 μs   (2.345 μs .. 2.546 μs)
-- std dev              336.4 ns   (248.6 ns .. 467.2 ns)
-- variance introduced by outliers: 93% (severely inflated)
--
-- Zipper; set node 150
-- benchmarking...
-- time                 8.838 μs   (8.768 μs .. 8.937 μs)
--                      0.999 R²   (0.998 R² .. 1.000 R²)
-- mean                 8.815 μs   (8.770 μs .. 8.895 μs)
-- std dev              205.4 ns   (116.7 ns .. 306.8 ns)
-- variance introduced by outliers: 25% (moderately inflated)
-- @
fnData :: String
fnData = "plh-multivariate.data"

main :: IO ()
main = do
  as <- getArgs

  case as of
    ["inspect"] -> do
      tr <- oneTree fn
      putStrLn $ "The tree has " <> show (length $ leaves tr) <> " leaves."
      let trRooted = either error id $ outgroup outgroups "root" tr
      L.putStrLn $ toNewick $ lengthToPhyloTree trRooted
      print $ getBranches trRooted
      let xs = itoList trRooted
      print xs
      let (pth, _) = fromMaybe (error "Gn_montanu not found.") $ ifind (\_ n -> n == "Gn_montanu") tr
      print pth
      let zp = zipper tr & downward (singular $ ix pth) & fromWithin _2 & focus .~ "BLAAAAAAAAAAAAAAAA" & rezip
      print zp
    ["read"] -> do
      putStrLn "Read trees."
      trs <- someTrees fn
      let l = length $ nub $ map T.fromLabeledTree trs
      unless (l == 1) (error "Trees have different topologies.")
      let trsRooted = map (either error id . outgroup outgroups "root") trs

      putStrLn "Get the posterior means and the posterior covariance matrix."
      let pm = getPosteriorMatrix trsRooted
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
