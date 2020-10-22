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

-- The source code formatter Ormolu messes up the comments describing the used
-- libraries.

{- ORMOLU_DISABLE -}

-- Global libraries.
import Control.Comonad
import Control.Lens hiding ((<.>))
import Control.Monad
import Criterion
import Data.Aeson
import Data.Bifunctor
import qualified Data.ByteString.Char8 as BS
import Data.Maybe
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import GHC.Generics
import Numeric.LinearAlgebra (Matrix, (<#), (<.>))
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import System.Environment
import System.Random.MWC hiding (uniform)

-- ELynx library providing functions to work with trees.
import ELynx.Tree

-- Mcmc library.
import Mcmc
import Mcmc.Tree

-- Local libraries provided together with this module.
import Calibration
import Constraint
import Tools
{- ORMOLU_ENABLE -}

-- File storing unrooted trees obtained from a Bayesian phylogenetic analysis.
-- The posterior means and covariances of the branch lengths are obtained from
-- these trees and used to approximate the phylogenetic likelihood.
fnInTrees :: FilePath
fnInTrees = "mcmc-examples/PhylogeneticLikelihoodMultivariate/data/plants_1.treelist.gz"

-- Base name of analysis.
bnAnalysis :: String
bnAnalysis = "plh-multivariate"

-- State space containing all parameters.
--
-- The topologies of the time and rate tree are equal. This is, however, not
-- ensured by the types. In the future, we may just use one tree storing both,
-- the times and the rates.
data I = I
  { -- Birth rate of time tree.
    _timeBirthRate :: Double,
    -- Death rate of time tree.
    _timeDeathRate :: Double,
    -- Height of the time tree in absolute time measured in million years, see
    -- calibrations.
    _timeHeight :: Double,
    -- Normalized time tree of height 1.0. Branch labels denote relative times;
    -- node labels denote relative node height.
    _timeTree :: Tree Double Double,
    -- Shape of the gamma distribution prior of the rates.
    _rateShape :: Double,
    -- Scale of the gamma distribution prior of the rates.
    _rateScale :: Double,
    -- Rate tree. Branch labels denote relative rates; node labels are unused.
    _rateTree :: Tree Double ()
    -- Remark: Let t and r be the lengths of a branch of the time and rate trees
    -- respectively. The length d of this branch measured in number of
    -- substitutions is d=t*r. Since the time tree is normalized, the time tree
    -- height is implicitly covered by r. The absolute rate is R = r/h, where h
    -- is the height of the tree. Similarly, absolute time is T = t*h.
    --
    -- I think this is a relatively clean solution. The absolute tree height is
    -- only determined by the calibrations, and not by the phylogenetic
    -- likelihood.
  }
  deriving (Generic)

-- Create accessors (lenses) to the parameters in the state space.Nothing
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
      _timeDeathRate = 1.0,
      _timeHeight = 1000.0,
      _timeTree = t',
      _rateShape = 10.0,
      _rateScale = 2.0,
      _rateTree = bimap (const 1.0) (const ()) t
    }
  where
    t' = extend rootHeight $ normalizeHeight $ makeUltrametric t

-- Calibrations are defined in the module 'Calibration'.

-- Constraints are defined in the module 'Constraint'.

-- Prior.
pr :: [Calibration] -> [Constraint] -> I -> Log Double
pr cb cs (I l m h t k th r) =
  product' $
    [ -- Exponential prior on the birth and death rates of the time tree.
      exponential 1 l,
      exponential 1 m,
      -- No prior on the height of the time tree but see the calibrations below.
      --
      -- Birth and death process prior of the time tree.
      birthDeath l m t,
      -- The reciprocal shape is exponentially distributed such that higher
      -- shape values are favored.
      exponential 10 k1,
      -- The scale is exponentially distributed.
      exponential 1 th,
      -- The prior of the branch-wise rates is gamma distributed.
      uncorrelatedGammaNoStem k th r
    ]
      ++ calibrations cb h t
      ++ constraints cs t
  where
    k1 = 1 / k

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
    times = getBranches (x ^. timeTree)
    rates = getBranches (x ^. rateTree)
    distances = sumFirstTwo $ V.zipWith (*) times rates

-- Proposals for the time tree.
proposalsTimeTree :: Show a => Tree e a -> [Proposal I]
proposalsTimeTree t =
  (timeTree @~ pulleyUltrametric 0.01 "Time tree root" 5 True) :
  [ (timeTree . nodeAt pth)
      @~ slideNodeUltrametric 0.01 ("Time tree node " ++ show lb) 1 True
    | (pth, lb) <- itoList t,
      -- Since the stem does not change the likelihood, it is set to zero, and
      -- we do not slide the root node.
      not (null pth),
      -- Also, we do not slide leaf nodes, since this would break
      -- ultrametricity.
      not (null $ forest $ current $ unsafeGoPath pth $ fromTree t)
  ]
    ++ [ (timeTree . nodeAt pth)
           @~ scaleSubTreeUltrametric 0.01 ("Time tree node " ++ show lb) 1 True
         | (pth, lb) <- itoList t,
           -- Don't scale the sub tree of the root node, because we are not
           -- interested in changing the length of the stem.
           not (null pth),
           -- Sub trees of leaves cannot be scaled.
           not (null $ forest $ current $ unsafeGoPath pth $ fromTree t)
       ]

-- Proposals for the rate tree.
proposalsRateTree :: Show a => Tree e a -> [Proposal I]
proposalsRateTree t =
  (rateTree @~ pulley 0.1 "Rate tree root" 5 True) :
  [ (rateTree . nodeAt pth)
      @~ slideBranch 0.1 ("Rate tree branch " ++ show lb) 1 True
    | (pth, lb) <- itoList t,
      -- Since the stem does not change the likelihood, it is set to zero, and
      -- we do not slide the stem.
      not (null pth)
  ]
    ++ [ (rateTree . nodeAt pth)
           @~ scaleTree 100 ("Rate tree node " ++ show lb) 1 True
         | (pth, lb) <- itoList t,
           -- Path does not lead to a leaf.
           not (null $ forest $ current $ unsafeGoPath pth $ fromTree t)
       ]

-- The complete cycle includes proposals for the other parameters.
ccl :: Show a => Tree e a -> Cycle I
ccl t =
  fromList $
    [ timeBirthRate @~ scaleUnbiased 10 "Time birth rate" 10 True,
      timeDeathRate @~ scaleUnbiased 10 "Time death rate" 10 True,
      timeHeight @~ scaleUnbiased 3000 "Time height" 10 True,
      rateShape @~ scaleUnbiased 10 "Rate shape" 10 True,
      rateScale @~ scaleUnbiased 10 "Rate scale" 10 True
    ]
      ++ proposalsTimeTree t
      ++ proposalsRateTree t

monParams :: [MonitorParameter I]
monParams =
  [ _timeBirthRate >$< monitorDouble "TimeBirthRate",
    _timeDeathRate >$< monitorDouble "TimeDeathRate",
    _timeHeight >$< monitorDouble "TimeHeight",
    _rateShape >$< monitorDouble "RateShape",
    _rateScale >$< monitorDouble "RateScale"
  ]

monStdOut :: MonitorStdOut I
monStdOut = monitorStdOut monParams 1

monFileParams :: MonitorFile I
monFileParams = monitorFile "-params" monParams 1

absoluteTimeTree :: I -> Tree Double Int
absoluteTimeTree s = identify $ first (* h) t
  where
    h = s ^. timeHeight
    t = s ^. timeTree

monFileTimeTree :: MonitorFile I
monFileTimeTree = monitorFile "-timetree" [absoluteTimeTree >$< monitorTree "TimeTree"] 1

monFileRateTree :: MonitorFile I
monFileRateTree = monitorFile "-ratetree" [_rateTree >$< monitorTree "RateTree"] 1

-- Collect monitors to standard output and files, as well as batch monitors.
mon :: Monitor I
mon = Monitor monStdOut [monFileParams, monFileTimeTree, monFileRateTree] []

-- Number of burn in iterations.
nBurnIn :: Maybe Int
-- nBurnIn = Just 30
nBurnIn = Just 3000

-- Auto tuning period.
nAutoTune :: Maybe Int
-- nAutoTune = Just 10
nAutoTune = Just 100

-- Number of Metropolis-Hasting iterations after burn in.
nIterations :: Int
-- nIterations = 10
nIterations = 10000

-- The rooted tree with posterior mean branch lengths will be stored in a file
-- with this name.
getMeanTree :: IO (Tree Double BS.ByteString)
getMeanTree = oneTree $ bnAnalysis ++ ".meantree"

-- Get the posterior branch length means, the inverted covariance matrix, and
-- the determinant of the covariance matrix.
getData :: IO (Vector Double, Matrix Double, Double)
getData = do
  (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' fnData
  let sigmaInv = L.fromRows sigmaInvRows
  return (mu, sigmaInv, logSigmaDet)
  where
    fnData = bnAnalysis ++ ".data"

main :: IO ()
main = do
  -- Get arguments.
  as <- getArgs
  case as of
    -- Read in all trees, calculate posterior means and covariances of the
    -- branch lengths, and find the midpoint root of the mean tree.
    ["prepare"] -> do
      prepare fnInTrees bnAnalysis
    -- Run the Metropolis-Hastings sampler.
    ["run"] -> do
      runMetropolisHastings
    -- Continue sampling.
    ["continue", n] -> do
      continueMetropolisHastings (read n)
    -- Print conversion of leaves.
    ["convert"] -> do
      convert
    -- Benchmark different functions used by the MCMC sampler.
    ["bench"] -> do
      runBench
    -- Inspect different objects; useful for debugging.
    ["inspect"] -> do
      inspect
    -- Print usage instructions if none of the previous commands was entered.
    _ -> putStrLn "Use one command of: [prepare|run|continue N|convert|bench|inspect]."

-- Run the Metropolis-Hastings sampler.
runMetropolisHastings :: IO ()
runMetropolisHastings = do
  -- Read the mean tree and the posterior means and covariances.
  meanTree <- getMeanTree
  (mu, sigmaInv, logSigmaDet) <- getData
  putStrLn "The posterior means of the branch lengths are:"
  print mu
  -- Initialize a starting state using the mean tree.
  let start = initWith $ identify meanTree
      pr' = pr (getCalibrations meanTree) (getConstraints meanTree)
      lh' = lh mu sigmaInv logSigmaDet
      ccl' = ccl $ identify meanTree
  -- Create a seed value for the random number generator. Actually, the
  -- 'create' function is deterministic, but useful during development. For
  -- real analyses, use 'createSystemRandom'.
  g <- create
  -- Construct the status of the Markov chain.
  let s =
        force . saveWith 1 . debug $
          -- Have a look at the 'status' function to understand the
          -- different parameters.
          status
            "plh-multivariate"
            pr'
            lh'
            ccl'
            mon
            start
            nBurnIn
            nAutoTune
            nIterations
            g
  -- Run the Markov chain.
  void $ mh s

continueMetropolisHastings :: Int -> IO ()
continueMetropolisHastings n = do
  meanTree <- getMeanTree
  (mu, sigmaInv, logSigmaDet) <- getData
  -- Load the MCMC status.
  s <-
    loadStatus
      (pr (getCalibrations meanTree) (getConstraints meanTree))
      (lh mu sigmaInv logSigmaDet)
      (ccl $ identify meanTree)
      mon
      "plh-multivariate.mcmc"
  void $ mhContinue n s

convert :: IO ()
convert = do
  meanTree <- getMeanTree
  let tbl = zip (map BS.unpack $ leaves meanTree) (map show $ leaves $ identify meanTree)
  putStr $ unlines $ map show tbl

-- Benchmark different functions used by the MCMC sampler.
runBench :: IO ()
runBench = do
  tr <- getMeanTree
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
      pr' = pr (getCalibrations tr) (getConstraints tr)
  benchmark $ nf pr' i
  (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
  let sigmaInv = L.fromRows sigmaInvRows
      lh' = lh mu sigmaInv logSigmaDet
  putStrLn "Benchmark calculation of likelihood."
  benchmark $ nf lh' i

  putStrLn "Benchmark identify."
  benchmark $ nf identify tr

-- Inspect different objects; useful for debugging.
inspect :: IO ()
inspect = do
  tr <- getMeanTree
  let lvs = leaves tr
  putStrLn $ "The mean tree has " <> show (length lvs) <> " leaves."

  let i = initWith $ identify tr
      pr' = pr (getCalibrations tr) (getConstraints tr)
  putStrLn $ "Test if time tree is ultrametric: " <> show (ultrametric $ i ^. timeTree)
  putStrLn $ "Initial prior: " <> show (pr' i) <> "."
  putStrLn "Load posterior means and covariances."
  (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
  let sigmaInv = L.fromRows sigmaInvRows
      lh' = lh mu sigmaInv logSigmaDet
  putStrLn $ "Initial log-likelihood: " <> show (ln $ lh' i) <> "."

  putStrLn "Paths to calibration nodes."
  putStrLn "Slight difference: 179, 183 younger than 245 247."
  putStrLn $ "Mrca of 179 and 183: " ++ show (mrca [179, 183] $ identify tr)
  putStrLn $ "Mrca of 245 and 247: " ++ show (mrca [245, 247] $ identify tr)
  putStrLn "More extreme difference: 34, 38 younger than 62, 64."
  putStrLn $ "Mrca of 34 and 38: " ++ show (mrca [34, 38] $ identify tr)
  putStrLn $ "Mrca of 62 and 64: " ++ show (mrca [62, 64] $ identify tr)

-- Post scriptum:
--
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
