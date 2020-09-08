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
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.ByteString.Char8 as BS
import Data.List
import Data.Maybe
import qualified Data.Set as S
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import qualified ELynx.Topology.Rooted as T
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
{- ORMOLU_ENABLE -}

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
    -- Time tree. Branch labels denote time; node labels denote node height.
    _timeTree :: Tree Double Double,
    -- Rate tree with mean branch length (1.0 / height of time tree). Branch
    -- labels denote rate; node labels are unused.
    _rateTree :: Tree Double ()
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
    { _timeBirthRate = hR,
      _timeDeathRate = hR,
      _timeTree = t',
      _rateTree = bimap (const hR) (const ()) t
    }
  where
    hI = 1000
    hR = recip hI
    t' = bimap (* hI) (* hI) $ extend rootHeight $ normalizeHeight $ makeUltrametric t

-- Calibrations are defined in 'Calibrations'. See 'calibratedNodes'.

-- Calibration prior with uniform bounds.
cals :: [Calibration] -> I -> [Log Double]
cals xs s = [calibrateUniformSoft 0.01 a b x t | (x, a, b) <- xs]
  where
    t = s ^. timeTree

-- Constraints define node orders.
--
-- (a, b) ensures node a to be younger than node b.
type Constraint = (Path, Path)

-- The constraint induced by a horizontal gene transfer does not contradict the
-- node order of the substitution-like tree obtained from the sequences.
constrainedNodes :: Tree e BS.ByteString -> [Constraint]
constrainedNodes t = [(young, old)]
  where
    young =
      fromMaybe
        (error "constrainedNodes: No MRCA young.")
        -- (mrca [213, 200])
        (mrca ["Pt_vittata", "Po_acrosti"] t)
    old =
      fromMaybe
        (error "constrainedNodes: No MRCA old.")
        -- (mrca [144, 143, 142])
        (mrca ["Me_tosanus", "Me_vincent", "No_aenigma"] t)

-- Constraint prior.
consts :: [Constraint] -> I -> [Log Double]
consts xs s =
  [constrainSoft 0.01 y o $ s ^. timeTree | (y, o) <- xs]

-- Prior.
pr :: [Calibration] -> [Constraint] -> I -> Log Double
pr cb cs s@(I l m t r) =
  product' $
    [ -- Exponential prior on the birth and death rates of the time tree.
      exponential 1 l,
      exponential 1 m,
      -- Birth and death process prior of the time tree.
      birthDeath l m t,
      -- The prior of the branch-wise rates is gamma distributed with mean
      -- (1.0/timeTreeHeight).
      branchesWith (gamma 1.0 (recip h)) r
    ]
      ++ cals cb s
      ++ consts cs s
  where
    h = label t

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
    times = getBranches (x ^. timeTree)
    rates = getBranches (x ^. rateTree)
    distances = sumFirstTwo $ V.zipWith (*) times rates

-- Slide node proposals for the time tree.
--
-- Since the stem does not change the likelihood, we do not slide the root node.
--
-- Also, we do not slide leaf nodes, since this would break ultrametricity.
proposalsTimeTree :: Show a => Tree e a -> [Proposal I]
proposalsTimeTree t =
  (timeTree @~ scaleTreeUltrametric "time tree scale" 10 3000 True) :
  [ (timeTree . nodeAt pth)
      @~ slideRootUltrametric ("time tree slide node " ++ show lb) 1 2.0 True
    | (pth, lb) <- itoList t,
      -- Path does not lead to the root.
      not (null pth),
      -- Path does not lead to a leaf.
      not (null $ forest $ current $ unsafeGoPath pth $ fromTree t)
  ]
    ++ [ (timeTree . nodeAt pth)
           @~ scaleSubTreeUltrametric ("time tree scale sub tree " ++ show lb) 1 2.0 True
         | (pth, lb) <- itoList t,
           -- Don't scale the sub tree of the root node, because we are not
           -- interested in the length of the stem.
           not (null pth),
           -- Sub trees of leaves cannot be scaled.
           not (null $ forest $ current $ unsafeGoPath pth $ fromTree t)
       ]

-- Slide branch proposals for the rate tree.
--
-- Since the stem does not change the likelihood, we do not slide the stem.
proposalsRateTree :: Show a => Tree e a -> [Proposal I]
proposalsRateTree t =
  [ (rateTree . nodeAt pth)
      @~ slideStem ("rate tree slide branch " ++ show lb) 1 0.001 True
    | (pth, lb) <- itoList t,
      -- Path does not lead to the root.
      not (null pth)
  ]
    ++ [ (rateTree . nodeAt pth)
           @~ scaleTree ("rate tree scale sub tree " ++ show lb) 1 1000 True
         | (pth, lb) <- itoList t,
           -- Path does not lead to a leaf.
           not (null $ forest $ current $ unsafeGoPath pth $ fromTree t)
       ]

-- Lens to the tree tuple; useful to create a contrary proposal scaling both
-- trees in opposite directions.
trLens :: Lens' I (Tree Double Double, Tree Double ())
trLens = lens (\x -> (_timeTree x, _rateTree x)) (\x (t, r) -> x {_timeTree = t, _rateTree = r})

-- The complete cycle includes slide proposals of higher weights for the other
-- parameters.
ccl :: Show a => Tree e a -> Cycle I
ccl t =
  fromList $
    [ timeBirthRate @~ scaleUnbiased "time birth rate" 10 10 True,
      timeDeathRate @~ scaleUnbiased "time death rate" 10 10 True,
      trLens @~ scaleTreesContrarily "time/rate tree contra scale" 10 3000 True
    ]
      ++ proposalsTimeTree t
      ++ proposalsRateTree t

monParams :: [MonitorParameter I]
monParams =
  [ _timeBirthRate @. monitorDouble "TimeBirthRate",
    _timeDeathRate @. monitorDouble "TimeDeathRate",
    (label . _timeTree) @. monitorDouble "TimeTreeRootHeight",
    view (timeTree . nodeAt [1,1] . rootLabel) @. monitorDouble "TimeTreeRootHeight"
  ]

monStdOut :: MonitorStdOut I
-- Do not monitor rateGammaShape to standard output because screen is not wide enough.
monStdOut = monitorStdOut monParams 1

monFileParams :: MonitorFile I
monFileParams = monitorFile "-params" monParams 1

monFileTimeTree :: MonitorFile I
monFileTimeTree = monitorFile "-timetree" [(identify . _timeTree) @. monitorTree "TimeTree"] 1

monFileRateTree :: MonitorFile I
monFileRateTree = monitorFile "-ratetree" [_rateTree @. monitorTree "RateTree"] 1

-- Collect monitors to standard output and files, as well as batch monitors.
mon :: Monitor I
mon = Monitor monStdOut [monFileParams, monFileTimeTree, monFileRateTree] []

-- Number of burn in iterations.
nBurnIn :: Maybe Int
nBurnIn = Just 3000
-- nBurnIn = Just 30

-- Auto tuning period.
nAutoTune :: Maybe Int
nAutoTune = Just 100
-- nAutoTune = Just 10

-- Number of Metropolis-Hasting iterations after burn in.
nIterations :: Int
nIterations = 10000
-- nIterations = 10

-- The posterior branch length means and covariances will be stored in a file
-- with this name.
fnData :: String
fnData = "plh-multivariate.data"

-- The rooted tree with posterior mean branch lengths will be stored in a file
-- with this name.
fnMeanTree :: FilePath
fnMeanTree = "plh-multivariate.meantree"

-- Read the mean tree and the posterior means and covariances.
readMeans :: IO (Tree Double BS.ByteString, Vector Double, Matrix Double, Double)
readMeans = do
  meanTree <- oneTree fnMeanTree
  (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
  let sigmaInv = L.fromRows sigmaInvRows
  return (meanTree, mu, sigmaInv, logSigmaDet)

main :: IO ()
main = do
  -- Get arguments.
  as <- getArgs

  case as of
    -- Read in all trees, calculate posterior means and covariances of the
    -- branch lengths, and find the midpoint root of the mean tree.
    ["read"] -> do
      putStrLn "Read trees; skip a burn in of 1000 trees."
      trs <- drop 1000 <$> someTrees fnTreeList

      -- Check if trees have the same topology.
      let l = length $ nub $ map T.fromLabeledTree trs
      unless (l == 1) (error "Trees have different topologies.")

      putStrLn "Calculate mean branch lengths."
      let pm = getPosteriorMatrix trs
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
      print $ toNewick $ measurableToPhyloTree tr
      putStrLn $ "Save the mean tree to " <> fnMeanTree <> "."
      BL.writeFile fnMeanTree (toNewick $ measurableToPhyloTree tr)

      putStrLn "Root the trees at the midpoint of the mean tree."
      let outgroups = fst $ fromBipartition $ either error id $ bipartition tr
          trsRooted = map (first fromLength . either error id . outgroup outgroups "root" . first Length) trs
      putStrLn "Get the posterior means and the posterior covariance matrix."
      let pmR = getPosteriorMatrixRooted trsRooted
          (mu, sigma) = L.meanCov pmR
          (sigmaInv, (logSigmaDet, _)) = L.invlndet $ L.unSym sigma

      putStrLn $ "Save the posterior means and covariances to " <> fnData <> "."
      encodeFile fnData (mu, L.toRows sigmaInv, logSigmaDet)

    -- Benchmark different functions used by the MCMC sampler.
    ["bench"] -> do
      tr <- oneTree fnMeanTree
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
          pr' = pr (calibratedNodes tr) (constrainedNodes tr)
      benchmark $ nf pr' i
      (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' "plh-multivariate.data"
      let sigmaInv = L.fromRows sigmaInvRows
          lh' = lh mu sigmaInv logSigmaDet
      putStrLn "Benchmark calculation of likelihood."
      benchmark $ nf lh' i

      putStrLn "Benchmark identify."
      benchmark $ nf identify tr

    -- Inspect different objects; useful for debugging.
    ["inspect"] -> do
      tr <- oneTree fnMeanTree
      let lvs = leaves tr
      putStrLn $ "The mean tree has " <> show (length lvs) <> " leaves."

      let i = initWith $ identify tr
          pr' = pr (calibratedNodes tr) (constrainedNodes tr)
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

    -- Run the Metropolis-Hastings sampler.
    ["run"] -> do
      -- Read the mean tree and the posterior means and covariances.
      (meanTree, mu, sigmaInv, logSigmaDet) <- readMeans
      putStrLn "The posterior means of the branch lengths are:"
      print mu
      -- Initialize a starting state using the mean tree.
      let start = initWith $ identify meanTree
          pr' = pr (calibratedNodes meanTree) (constrainedNodes meanTree)
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
    ["continue", n] -> do
      (meanTree, mu, sigmaInv, logSigmaDet) <- readMeans
      -- Load the MCMC status.
      s <-
        loadStatus
          (pr (calibratedNodes meanTree) (constrainedNodes meanTree))
          (lh mu sigmaInv logSigmaDet)
          (ccl $ identify meanTree)
          mon
          "plh-multivariate.mcmc"
      void $ mhContinue (read n) s
    -- Print usage instructions if none of the previous commands was entered.
    ["convert"] -> do
      -- Print conversion of leaves.
      (meanTree, _, _, _) <- readMeans
      let tbl = zip (map BS.unpack $ leaves meanTree) (map show $ leaves $ identify meanTree)
      putStr $ unlines $ map show tbl
    _ -> putStrLn "Use one command of: [read|bench|inspect|run|continue N|convert]!"

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
