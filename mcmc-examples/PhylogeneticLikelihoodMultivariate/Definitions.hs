{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Definitions
-- Description :  State space, prior function, likelihood function and more.
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Oct 22 20:05:05 2020.
module Definitions
  ( fnInTrees,
    bnAnalysis,
    I (..),
    initWith,
    priorFunction,
    likelihoodFunction,
    proposals,
    monitor,
    burnInSpec,
    nIterations,
  )
where

import Control.Lens
import Data.Aeson
import Data.Bifunctor
import qualified Data.Vector.Storable as V
import GHC.Generics
import qualified Numeric.LinearAlgebra as L
import Numeric.Log

-- import Debug.Trace

-- Disable the syntax formatter Ormolu to highlight relevant module imports.
{- ORMOLU_DISABLE -}
-- The ELynx library includes functions to work on trees.
import ELynx.Tree

-- The Mcmc library includes the Metropolis-Hastings sampler.
import Mcmc
import Mcmc.Tree

-- Local modules.
import Calibrations
import Constraints
import Tools
{- ORMOLU_ENABLE -}

-- | File storing unrooted trees obtained from a Bayesian phylogenetic analysis.
-- The posterior means and covariances of the branch lengths are obtained from
-- these trees and used to approximate the phylogenetic likelihood.
fnInTrees :: FilePath
fnInTrees = "mcmc-examples/PhylogeneticLikelihoodMultivariate/data-landplants/alignment_201020.phy_gtr_1.treelist"

-- | Base name of analysis.
bnAnalysis :: String
bnAnalysis = "plh-multivariate"

-- | State space containing all parameters.
--
-- We are interested in inferring an ultrametric tree with branch lengths
-- measured in units of time (e.g., in million years). Let T be the length of a
-- branch of the time tree, and R be the absolute evolutionary rate on this
-- branch. Then, the length of this very same branch measured in average number
-- of substitutions is d=T*R.
--
-- Internally, a relative time t and relative rate r are used such that the
-- branch length measured in average number of substitution is
-- d=T*R=(t*h)*(r*mu), where h is the root height of the time tree, and mu is
-- the mean rate.
--
-- In brief, the relative time and rate are defined as t=T/h, and r=R/mu, where
-- h is the root height and mu is the mean rate.
--
-- This has two advantages:
--
-- 1. The ultrametric tree object storing the relative times is a normalized
--    tree with root height 1.0.
--
-- 2. The relative rates have a mean of 1.0.
--
-- Remark: The topologies of the time and rate tree are equal. This is, however,
-- not ensured by the types. One could just use one tree storing both, the times
-- and the rates.
data I = I
  { -- | Birth rate of relative time tree.
    _timeBirthRate :: Double,
    -- | Death rate of relative time tree.
    _timeDeathRate :: Double,
    -- | Height of the absolute time tree in unit time. Here, we use units of
    -- million years; see the calibrations.
    _timeHeight :: Double,
    -- | Normalized time tree of height 1.0. Branch labels denote relative
    -- times; node labels store relative node height and names.
    _timeTree :: Tree Length (HeightLabel Name),
    -- | The mean of the absolute rates.
    _rateMean :: Double,
    -- | The shape of the relative rates.
    _rateShape :: Double,
    -- | Relative rate tree. Branch labels denote relative rates with mean 1.0;
    -- node labels store names.
    _rateTree :: Tree Length Name
  }
  deriving (Generic)

-- Create accessors (lenses) to the parameters in the state space.
makeLenses ''I

-- Allow storage of the trace as JSON.
instance ToJSON I

instance FromJSON I

-- | Initial state.
--
-- The given tree is used to initiate the time and rate trees. For the time
-- tree, the terminal branches are elongated such that the tree becomes
-- ultrametric ('makeUltrametric'). For the rate tree, we just use the topology
-- and set all rates to 1.0.
initWith :: Tree Length Name -> I
initWith t =
  I
    { _timeBirthRate = 1.0,
      _timeDeathRate = 1.0,
      _timeHeight = 1000.0,
      _timeTree = t',
      _rateMean = 1 / 1000.0,
      _rateShape = 10,
      _rateTree = setStem 0 $ first (const 1.0) t
    }
  where
    t' = toHeightTree $ normalizeHeight $ makeUltrametric t

-- The calibrations are defined in the module 'Calibration'.

-- The constraints are defined in the module 'Constraint'.

-- | Prior function.
priorFunction :: [Calibration] -> [Constraint] -> PriorFunction I
priorFunction cb cs (I l m h t mu k r) =
  product' $
    [ -- Exponential prior on the birth and death rates of the time tree.
      exponential 1 l,
      exponential 1 m,
      -- No explicit prior on the height of the time tree. However, calibrations
      -- are used (see below).
      --
      -- Birth and death process prior on the time tree.
      birthDeath WithoutStem l m 1.0 t,
      -- Exponential prior on the rate mean.
      exponential 1 mu,
      -- Exponential prior on the shape parameter of the rate prior.
      exponential 1 k,
      -- Uncorrelated log normal prior on the branch-wise rates.
      uncorrelatedGamma WithoutStem k k1 r
    ]
      -- Add the calibrations and constraint.
      ++ calibrations cb h t
      ++ constraints cs t
  where
    k1 = recip k

-- Log of density of multivariate normal distribution with given parameters.
-- https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
--
-- The constant @k * log (2*pi)@ was left out on purpose.
logDensityMultivariateNormal ::
  -- Mean vector.
  V.Vector Double ->
  -- Inverted covariance matrix.
  L.Matrix Double ->
  -- Log of determinant of covariance matrix.
  Double ->
  -- Value vector.
  V.Vector Double ->
  Log Double
logDensityMultivariateNormal mu sigmaInv logSigmaDet xs =
  Exp $ (-0.5) * (logSigmaDet + ((dxs L.<# sigmaInv) L.<.> dxs))
  where
    dxs = xs - mu

-- | Approximation of the phylogenetic likelihood using a multivariate normal
-- distribution.
likelihoodFunction ::
  -- | Mean vector.
  V.Vector Double ->
  -- | Inverted covariance matrix.
  L.Matrix Double ->
  -- | Log of determinant of covariance matrix.
  Double ->
  LikelihoodFunction I
likelihoodFunction mu sigmaInv logSigmaDet x =
  logDensityMultivariateNormal mu sigmaInv logSigmaDet distances
  where
    times = getBranches (x ^. timeTree)
    rates = getBranches (x ^. rateTree)
    tH = x ^. timeHeight
    rMu = x ^. rateMean
    distances = V.map (* (tH * rMu)) $ sumFirstTwo $ V.zipWith (*) times rates

-- Proposals for the time tree.
proposalsTimeTree :: Show a => Tree e a -> [Proposal I]
proposalsTimeTree t =
  -- Pulley on the root node.
  (timeTree @~ pulleyUltrametric t 0.1 (PName "Time tree root") (PWeight 10) Tune) :
  -- Slide nodes excluding the root and the leaves.
  [ {-# SCC slideNodeUltrametric #-}
    (timeTree . subTreeAtE pth)
      @~ slideNodeUltrametric 0.1 (PName $ "Time tree node " ++ show lb) (PWeight 1) Tune
    | (pth, lb) <- itoList $ identify t,
      -- Since the stem does not change the likelihood, it is set to zero, and
      -- we do not slide the root node.
      not (null pth),
      -- Also, we do not slide leaf nodes, since this would break
      -- ultrametricity.
      not $ null $ t ^. subTreeAtE pth . forestL
  ]
    -- Scale sub trees of inner nodes excluding the root and the leaves.
    ++ [ {-# SCC scaleSubTreeUltrametric #-}
         (timeTree . subTreeAtE pth)
           @~ scaleSubTreeUltrametric s 0.1 (PName $ "Time tree node " ++ show lb) (PWeight 1) Tune
         | (pth, lb) <- itoList $ identify t,
           let s = t ^. subTreeAtE pth,
           -- Don't scale the sub tree of the root node, because we are not
           -- interested in changing the length of the stem.
           not $ null pth,
           -- Sub trees of leaves cannot be scaled.
           not $ null $ forest s
       ]

-- Proposals for the rate tree.
proposalsRateTree :: Show a => Tree e a -> [Proposal I]
proposalsRateTree t =
  -- Pulley on the root node.
  (rateTree @~ pulley 0.1 (PName "Rate tree root") (PWeight 10) Tune) :
  -- Scale branches excluding the stem.
  [ {-# SCC slideBranch #-}
    (rateTree . subTreeAtE pth)
      @~ scaleBranch 0.1 (PName $ "Rate tree branch " ++ show lb) (PWeight 1) Tune
    | (pth, lb) <- itoList $ identify t,
      -- Since the stem does not change the likelihood, it is set to zero, and
      -- we do not slide the stem.
      not (null pth)
  ]
    -- Scale trees of inner nodes excluding the root and the leaves.
    ++ [ {-# SCC scaleTree #-}
         (rateTree . subTreeAtE pth)
           @~ scaleTree s WithoutStem 100 (PName $ "Rate tree node " ++ show lb) (PWeight 1) Tune
         | (pth, lb) <- itoList $ identify t,
           let s = t ^. subTreeAtE pth,
           -- Path does not lead to a leaf.
           not $ null $ forest s
       ]

-- Lens for a contrary proposal.
timeHeightRateNormPair :: Lens' I (Double, Double)
timeHeightRateNormPair =
  lens
    (\x -> (x ^. timeHeight, x ^. rateMean))
    (\x (h, mu) -> x {_timeHeight = h, _rateMean = mu})

-- | The proposal cycle includes proposals for the other parameters.
proposals :: Show a => Tree e a -> Cycle I
proposals t =
  fromList $
    [ timeBirthRate @~ scaleUnbiased 10 (PName "Time birth rate") (PWeight 20) Tune,
      timeDeathRate @~ scaleUnbiased 10 (PName "Time death rate") (PWeight 20) Tune,
      timeHeight @~ scaleUnbiased 3000 (PName "Time height") (PWeight 20) Tune,
      rateMean @~ scaleUnbiased 10 (PName "Rate mean") (PWeight 20) Tune,
      rateShape @~ scaleUnbiased 10 (PName "Rate shape") (PWeight 20) Tune,
      timeHeightRateNormPair @~ scaleContrarily 10 0.1 (PName "Time height, rate mean") (PWeight 20) Tune
    ]
      ++ proposalsTimeTree t
      ++ proposalsRateTree t

-- -- Monitor the average rate. Useful, because it should not deviate from 1.0 too
-- -- much.
-- getAverageBranchLength :: (I -> Tree Length a) -> I -> Double
-- getAverageBranchLength f x = (/ n) $ fromLength $ totalLength r
--   where
--     r = f x
--     n = fromIntegral $ length r - 1

-- -- 0 if the height correctly calculated?
-- monDeltaHeight :: Path -> I -> Double
-- monDeltaHeight pth x = fromLength (t ^. labelL . heightL - rootHeight t)
--   where t = x ^. timeTree . subTreeAtE pth

-- Monitor parameters.
monParams :: [MonitorParameter I]
monParams =
  [ _timeBirthRate >$< monitorDouble "TimeBirthRate",
    _timeDeathRate >$< monitorDouble "TimeDeathRate",
    _timeHeight >$< monitorDouble "TimeHeight",
    _rateMean >$< monitorDouble "RateMean",
    _rateShape >$< monitorDouble "RateShape"
    -- getAverageBranchLength _rateTree >$< monitorDouble "RateAverage",
    -- monDeltaHeight [0,0] >$< monitorDoubleE "Delta Root"
  ]

-- Monitor to standard output.
monStdOut :: MonitorStdOut I
-- Screen width doesn't support more than four parameter monitors.
monStdOut = monitorStdOut (take 4 monParams) 1

-- Get the height of the node at path. Useful to have a look at calibrated nodes.
getTimeTreeNodeHeight :: Path -> I -> Double
getTimeTreeNodeHeight p x = (* h) $ fromLength $ t ^. subTreeAtE p . labelL . heightL
  where
    t = x ^. timeTree
    h = x ^. timeHeight

-- Monitor the height of calibrated nodes.
monCalibratedNodes :: [Calibration] -> [MonitorParameter I]
monCalibratedNodes cb = [getTimeTreeNodeHeight p >$< monitorDouble (nm n a b) | (n, p, a, b) <- cb]
  where
    nm s l r = "Calibration " ++ s ++ " (" ++ show l ++ ", " ++ show r ++ ")"

-- Get the difference in height of the nodes at path. Useful to have a look at
-- constrained nodes. Positive if constraint is honored.
getTimeTreeDeltaNodeHeight :: Path -> Path -> I -> Double
getTimeTreeDeltaNodeHeight y o x = getTimeTreeNodeHeight o x - getTimeTreeNodeHeight y x

-- Monitor the heights of constrained nodes.
monConstrainedNodes :: [Constraint] -> [MonitorParameter I]
monConstrainedNodes cs = [getTimeTreeDeltaNodeHeight y o >$< monitorDouble (nm n) | (n, y, o) <- cs]
  where
    nm s = "Constraint " ++ s

-- The file monitor is more verbose.
monFileParams :: [Calibration] -> [Constraint] -> MonitorFile I
monFileParams cb cs =
  monitorFile
    "-params"
    ( monParams
        ++ monCalibratedNodes cb
        ++ monConstrainedNodes cs
    )
    1

-- Monitor the time tree with absolute branch lengths, because they are more
-- informative.
absoluteTimeTree :: I -> Tree Length Name
absoluteTimeTree s = first (* h) t
  where
    h = either error id $ toLength $ s ^. timeHeight
    t = fromHeightTree $ s ^. timeTree

-- The time tree with absolute branch lengths is written to a separate file.
monFileTimeTree :: MonitorFile I
monFileTimeTree = monitorFile "-timetree" [absoluteTimeTree >$< monitorTree "TimeTree"] 1

-- The rate tree with relative rates is written to a separate file.
monFileRateTree :: MonitorFile I
monFileRateTree = monitorFile "-ratetree" [_rateTree >$< monitorTree "RateTree"] 1

-- | Monitor to standard output and files. Do not use any batch monitors for now.
monitor :: [Calibration] -> [Constraint] -> Monitor I
monitor cb cs = Monitor monStdOut [monFileParams cb cs, monFileTimeTree, monFileRateTree] []

-- | Number of burn in iterations and auto tuning period.
burnInSpec :: BurnIn
burnInSpec = BurnInWithAutoTuning 3000 100

-- | Number of Metropolis-Hasting iterations after burn in.
nIterations :: Int
nIterations = 40000
