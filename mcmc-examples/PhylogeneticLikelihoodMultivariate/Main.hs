{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

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
import qualified Data.ByteString.Conversion.From as L
import qualified Data.ByteString.Lazy.Char8 as L
import Data.List
import Data.Maybe
import qualified Data.Text.Lazy.Builder as B
import qualified Data.Set as S
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import qualified ELynx.Data.Topology.Rooted as T
import ELynx.Data.Tree
import ELynx.Export.Tree.Newick
import GHC.Generics
import Mcmc
import Numeric.LinearAlgebra (Matrix, (<#), (<.>))
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import Prior
import Proposal
import System.Environment
import System.Random.MWC
import Tree

-- State space containing all parameters.
--
-- The topologies of the time and rate tree are equal. This is, however, not
-- ensured by the types. In the future, we may just use one tree storing both,
-- the times and the rates.
data I = I
  { -- Birth rate parameter of time tree.
    _timeBirthRate :: Double,
    -- Height of root node measured in units of time.
    _timeRootHeight :: Double,
    -- Shape parameter k of gamma distribution of rate parameters. The scale
    -- parameter is determined such that the mean of the gamma distribution is
    -- 1.
    _rateGammaShape :: Double,
    -- Global mean of rate parameters.
    _rateMean :: Double,
    -- Time tree.
    _timeTree :: Tree Double Int,
    -- Rate tree.
    _rateTree :: Tree Double Int
  }
  deriving (Generic)

-- Create accessors (lenses) to the parameters in the state space.
makeLenses ''I

-- Allow storage of the trace as JSON.
instance ToJSON I

instance FromJSON I

-- Initial state.
--
-- The given tree is used to initiate the time and rate trees. For the time
-- tree, the terminal branches are elongated such that the tree becomes
-- ultrametric ('makeUltrametric'). For the rate tree, we just use the topology
-- and set all rates to 1.0.
initWith :: Tree Double Int -> I
initWith t =
  I
    { _timeBirthRate = 1.0,
      _timeRootHeight = height t',
      _rateGammaShape = 1.0,
      _rateMean = 1.0,
      _timeTree = t',
      _rateTree = first (const 1.0) t
    }
  where
    t' = makeUltrametric t

-- Prior.
pr :: I -> Log Double
pr (I l h k m t r) =
  product'
    [ -- Exponential prior on the birth rate of the time tree.
      exponentialWith 1.0 l,
      -- Exponential prior on the root height of the time tree.
      exponentialWith 10.0 h,
      -- Exponential prior on the shape of the gamma distribution.
      exponentialWith 10.0 k,
      -- The rate mean prior is a gamma distribution with shape k and mean one.
      gammaWith k (1 / k) m,
      -- Birth process prior on the branches of the time tree.
      branchesWith (exponentialWith l) t,
      -- The prior of the branch-wise rates is also gamma distributed.
      branchesWith (gammaWith k (1 / k)) r
    ]

-- File storing unrooted trees obtained from a Bayesian phylogenetic analysis.
-- The posterior means and covariances of the branch lengths are obtained from
-- these trees and used to approximate the phylogenetic likelihood.
fnTreeList :: FilePath
fnTreeList = "mcmc-examples/PhylogeneticLikelihoodMultivariate/data/plants_1.treelist.gz"

-- Helper function. Get all branches of a rooted tree. Store the branches in a
-- vector such that the two branches leading to the root are the first two
-- entries of the vector. Ignore the root branch.
getBranches :: Tree Double a -> Vector Double
getBranches (Node _ _ [l, r]) = V.fromList $ head ls : head rs : tail ls ++ tail rs
  where
    ls = branches l
    rs = branches r
getBranches _ = error "getBranches: Root node is not bifurcating."

-- Sum the first two elements of a vector. Needed to merge the two branches
-- leading to the root.
sumFirstTwo :: Vector Double -> Vector Double
sumFirstTwo v = (v V.! 0 + v V.! 1) `V.cons` V.drop 2 v

-- Get the posterior matrix of branch lengths of rooted trees. Merge the two
-- branch lengths leading to the root.
getPosteriorMatrixRooted :: [Tree Double a] -> Matrix Double
getPosteriorMatrixRooted = L.fromRows . map (sumFirstTwo . getBranches)

-- Get the posterior matrix of branch lengths of unrooted trees (trees with
-- multifurcating root nodes). Before midpoint rooting, the mean branch lengths
-- of the unrooted trees have to be determined.
getPosteriorMatrix :: [Tree Double a] -> Matrix Double
getPosteriorMatrix = L.fromRows . map (V.fromList . branches)

-- Log of density of multivariate normal distribution with given parameters.
-- https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
--
-- The constant @k * log (2*pi)@ was left out on purpose.
logDensityMultivariateNormal ::
  -- Mean vector.
  Vector Double ->
  -- Inverted covariance matrix.
  Matrix Double ->
  -- Log of determinant of covariance matrix.
  Double ->
  -- Value vector.
  Vector Double ->
  Log Double
logDensityMultivariateNormal mu sigmaInv logSigmaDet xs =
  Exp $ (-0.5) * (logSigmaDet + ((dxs <# sigmaInv) <.> dxs))
  where
    dxs = xs - mu

-- Phylogenetic likelihood using a multivariate normal distribution.
lh ::
  -- Mean vector.
  Vector Double ->
  -- Inverted covariance matrix.
  Matrix Double ->
  -- Log of determinant of covariance matrix.
  Double ->
  -- Current state.
  I ->
  Log Double
lh mu sigmaInv logSigmaDet x = logDensityMultivariateNormal mu sigmaInv logSigmaDet distances
  where
    times = getBranches $ x ^. timeTree
    rates = getBranches $ x ^. rateTree
    multiplier = x ^. timeRootHeight * x ^. rateMean
    distances = sumFirstTwo $ V.map (* multiplier) $ V.zipWith (*) times rates

-- Slide node proposals for the time tree.
--
-- Since the stem does not change the likelihood, we do not slide the root node.
--
-- Also, we do not slide leaf nodes, since this would break ultrametricity.
proposalsTimeTree :: Show a => Tree e a -> [Proposal I]
proposalsTimeTree t =
  [ timeTree >>> slideNode pth (n lb) 1
    | (pth, lb) <- itoList t,
      -- Path does not lead to the root.
      not (null pth),
      -- Path does not lead to a leaf.
      not (null $ forest $ current $ unsafeGoPath pth $ fromTree t)
  ]
  where
    n x = "timetree slide node " ++ show x

-- Slide branch proposals for the rate tree.
--
-- Since the stem does not change the likelihood, we do not slide the stem.
proposalsRateTree :: Show a => Tree e a -> [Proposal I]
proposalsRateTree t =
  [ rateTree >>> slideBranch pth (n lb) 1 1.0 True
    | (pth, lb) <- itoList t,
      -- Path does not lead to the root.
      not (null pth)
  ]
  where
    n x = "ratetree slide branch " ++ show x

-- The complete cycle includes slide proposals of higher weights for the other
-- parameters.
ccl :: Show a => Tree e a -> Cycle I
ccl t =
  fromList $
    [ timeBirthRate >>> slideSymmetric "time birth rate" 4 0.5 True,
      timeRootHeight >>> slideSymmetric "time root height" 4 0.5 True,
      rateGammaShape >>> slideSymmetric "rate gamma shape" 4 0.5 True,
      rateMean >>> slideSymmetric "rate mean" 4 0.5 True
    ]
      ++ proposalsTimeTree t
      ++ proposalsRateTree t

-- TODO: Provide more elaborate monitors.

-- TODO: Move this function into the library.

-- For now, only monitor the time tree.
monTimeTree :: MonitorParameter I
monTimeTree =
  MonitorParameter
    "timeTree"
    ( \x ->
        B.fromText $
          f $
            toNewick $
              lengthToPhyloTree $ first Length $ x ^. timeTree
    )
  where
    f s = fromMaybe (error "conversion failed") $ L.fromByteString $ L.toStrict s

-- Collect monitors to standard output and files, as well as batch monitors.
mon :: Monitor I
mon = Monitor (monitorStdOut [monTimeTree] 1) [] []

-- Number of burn in iterations.
nBurnIn :: Maybe Int
nBurnIn = Just 800

-- Auto tuning period.
nAutoTune :: Maybe Int
nAutoTune = Just 200

-- Number of Metropolis-Hasting iterations after burn in.
nIterations :: Int
nIterations = 1000

-- The posterior branch length means and covariances will be stored in a file
-- with this name.
fnData :: String
fnData = "plh-multivariate.data"

-- The rooted tree with posterior mean branch lengths will be stored in a file
-- with this name.
fnMeanTree :: FilePath
fnMeanTree = "plh-multivariate.meantree"

main :: IO ()
main = do
  -- Get arguments.
  as <- getArgs
  case as of
    ["read"] -> do
      putStrLn "Read trees; skip a burn in of 1000 trees."
      trs <- drop 1000 <$> someTrees fnTreeList
      let l = length $ nub $ map T.fromLabeledTree trs
      unless (l == 1) (error "Trees have different topologies.")

      putStrLn "Calculate mean branch lengths."
      let pm = getPosteriorMatrix $ map (first fromLength) trs
          (means, _) = L.meanCov pm
      putStrLn "The mean branch lengths are:"
      print means
      let meanTree =
            fromMaybe (error "Could not label tree with mean branch lengths") $
              setBranches (map Length $ V.toList means) (head trs)
          lvs = leaves meanTree
          trOutgroup = either error id $ outgroup (S.singleton $ head lvs) "root" meanTree
          tr = either error id $ midpoint trOutgroup
      putStrLn "The tree with mean branch lengths rooted at the midpoint:"
      print $ toNewick $ lengthToPhyloTree tr
      putStrLn $ "Save the mean tree to " <> fnMeanTree <> "."
      L.writeFile fnMeanTree (toNewick $ lengthToPhyloTree tr)

      putStrLn "Root the trees at the midpoint of the mean tree."
      let outgroups = fst $ fromBipartition $ either error id $ bipartition tr
          trsRooted = map (either error id . outgroup outgroups "root") trs

      putStrLn "Get the posterior means and the posterior covariance matrix."
      let pmR = getPosteriorMatrixRooted $ map (first fromLength) trsRooted
          (mu, sigma) = L.meanCov pmR
          (sigmaInv, (logSigmaDet, _)) = L.invlndet $ L.unSym sigma

      putStrLn $ "Save the posterior means and covariances to " <> fnData <> "."
      encodeFile fnData (mu, L.toRows sigmaInv, logSigmaDet)
    ["bench"] -> do
      tr <- first fromLength <$> oneTree fnMeanTree
      let (pth, _) =
            fromMaybe (error "Gn_montanu not found.") $
              ifind (\_ n -> n == "Gn_montanu") tr
      putStrLn $ "The path to \"Gn_montanu\" is: " <> show pth
      let bf1 =
            toTree . insertLabel "Bla"
              . fromMaybe (error "Path does not lead to a leaf.")
              . goPath pth
              . fromTree
      putStrLn $ "Change a leaf: " <> show (bf1 tr) <> "."
      putStrLn "Benchmark change a leaf."
      benchmark $ nf bf1 tr
      let bf2 =
            label . current
              . fromMaybe (error "Path does not lead to a leaf.")
              . goPath pth
              . fromTree
      putStrLn $ "Leaf to get: " <> show (bf2 tr) <> "."
      putStrLn "Benchmark get a leaf."
      benchmark $ nf bf2 tr
      putStrLn "Benchmark calculation of prior."
      let i = initWith $ identify tr
      benchmark $ nf pr i
      (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
      let sigmaInv = L.fromRows sigmaInvRows
          lh' = lh mu sigmaInv logSigmaDet
      putStrLn "Benchmark calculation of likelihood."
      benchmark $ nf lh' i
    ["inspect"] -> do
      tr <- first fromLength <$> oneTree fnMeanTree
      let lvs = leaves tr
      putStrLn $ "The tree has " <> show (length lvs) <> " leaves."
      let i = initWith $ identify tr
      putStrLn $ "Test if time tree is ultrametric: " <> show (ultrametric $ i ^. timeTree)
      putStrLn $ "Initial prior: " <> show (pr i) <> "."
      putStrLn "Load posterior means and covariances."
      (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
      let sigmaInv = L.fromRows sigmaInvRows
          lh' = lh mu sigmaInv logSigmaDet
      putStrLn $ "Initial log-likelihood: " <> show (ln $ lh' i) <> "."
    -- let tr' = fromMaybe (error "Could not drop leaves.") $ dropLeavesOnly tr
    -- putStrLn "Tree without leaves:"
    -- L.putStrLn $ toNewick $ lengthToPhyloTree $ bimap Length (L.pack . show) $ identify tr'
    -- putStrLn $ "Number of leaves: " ++ show (length $ leaves tr')
    ["run"] -> do
      tr <- first fromLength <$> oneTree fnMeanTree
      (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
      let sigmaInv = L.fromRows sigmaInvRows
      putStrLn "The posterior means of the branch lengths are:"
      print mu
      let start = initWith $ identify tr
      g <- create
      putStrLn "Construct status of the chain."
      let s =
            force $
              status
                "plh-multivariate"
                pr
                (lh mu sigmaInv logSigmaDet)
                (ccl tr)
                mon
                start
                nBurnIn
                nAutoTune
                nIterations
                g
      void $ mh s
    _ -> putStrLn "read|bench|inspect|run"

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
