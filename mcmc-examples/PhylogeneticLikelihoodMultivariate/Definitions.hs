{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Definitions
-- Description :  State space, prior, likelihood etc.
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
    cleaner,
    initWith,
    priorDistribution,
    likelihoodFunction,
    proposals,
    monitor,
    nBurnIn,
    nAutoTune,
    nIterations,
  )
where

import Control.Lens
import Data.Aeson
import Data.Bifunctor
import qualified Data.ByteString.Char8 as BS
import qualified Data.Vector.Storable as V
import GHC.Generics
import qualified Numeric.LinearAlgebra as L
import Numeric.Log

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

-- import Debug.Trace

-- | File storing unrooted trees obtained from a Bayesian phylogenetic analysis.
-- The posterior means and covariances of the branch lengths are obtained from
-- these trees and used to approximate the phylogenetic likelihood.
fnInTrees :: FilePath
fnInTrees = "mcmc-examples/PhylogeneticLikelihoodMultivariate/data/alignment_201020.phy_gtr_1.treelist"

-- | Base name of analysis.
bnAnalysis :: String
bnAnalysis = "plh-multivariate"

-- | State space containing all parameters.
--
-- We are interested in inferring an ultrametric tree with branch lengths
-- measured in unit time (e.g., in million years). Let T be the length of a
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
    _timeTree :: TimeTree,
    -- | The mean of the absolute rates.
    _rateMean :: Double,
    -- | The variance of the logarithm of the relative rates.
    _rateVariance :: Double,
    -- | Relative rate tree. Branch labels denote relative rates with mean 1.0;
    -- node labels store names.
    _rateTree :: RateTree
  }
  deriving (Generic)

-- Create accessors (lenses) to the parameters in the state space.
makeLenses ''I

-- Allow storage of the trace as JSON.
instance ToJSON I

instance FromJSON I

-- See 'cleaner'. This function makes the tree ultrametric again, normalizes the
-- tree and sets the height values accordingly.
cleanTimeTree :: I -> I
cleanTimeTree = timeTree %~ (toTimeTree . normalizeHeight . makeUltrametric . fromTimeTree)

-- | Clean the state periodically. Otherwise, the tree diverges from being
-- ultrametric.
cleaner :: Cleaner I
cleaner = Cleaner 100 cleanTimeTree

-- | Initial state.
--
-- The given tree is used to initiate the time and rate trees. For the time
-- tree, the terminal branches are elongated such that the tree becomes
-- ultrametric ('makeUltrametric'). For the rate tree, we just use the topology
-- and set all rates to 1.0.
initWith :: SubstitutionTree -> I
initWith t =
  I
    { _timeBirthRate = 1.0,
      _timeDeathRate = 1.0,
      _timeHeight = 1200.0,
      _timeTree = t',
      _rateMean = 1 / 1200.0,
      _rateVariance = 4,
      _rateTree = first (const 1.0) t
    }
  where
    t' = toTimeTree $ normalizeHeight $ makeUltrametric t

-- Calibrations are defined in the module 'Calibration'.

-- Constraints are defined in the module 'Constraint'.

-- | Prior distribution.
priorDistribution :: [Calibration] -> [Constraint] -> I -> Log Double
priorDistribution cb cs (I l m h t mu v r) =
  product' $
    [ -- Exponential prior on the birth and death rates of the time tree.
      exponential 1 l,
      exponential 1 m,
      -- No explicit prior on the height of the time tree. However, calibrations
      -- are used (see below).
      --
      -- Birth and death process prior on the time tree.
      birthDeath l m t,
      -- Gamma prior on the rate mean. The mean is k*theta=1/1200.
      gamma (recip $ 1200 * 5) 5 mu,
      -- Strong gamma prior on the rate variance. The mean is k*theta=0.01.
      gamma 1 0.01 v,
      -- Uncorrelated log normal prior on the branch-wise rates.
      uncorrelatedLogNormalNoStem 1.0 v r
    ]
      ++ calibrations cb h t
      ++ constraints cs t

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
  -- | Current state.
  I ->
  Log Double
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
  (timeTree @~ pulleyUltrametric 0.01 "Time tree root" 5 True) :
  [ (timeTree . subTreeAt pth)
      @~ slideNodeUltrametric 0.01 ("Time tree node " ++ show lb) 1 True
    | (pth, lb) <- itoList $ identify t,
      -- Since the stem does not change the likelihood, it is set to zero, and
      -- we do not slide the root node.
      not (null pth),
      -- Also, we do not slide leaf nodes, since this would break
      -- ultrametricity.
      not (null $ forest $ getSubTreeUnsafe pth t)
  ]
    ++ [ (timeTree . subTreeAt pth)
           @~ scaleSubTreeUltrametric 0.01 ("Time tree node " ++ show lb) 1 True
         | (pth, lb) <- itoList $ identify t,
           -- Don't scale the sub tree of the root node, because we are not
           -- interested in changing the length of the stem.
           not (null pth),
           -- Sub trees of leaves cannot be scaled.
           not (null $ forest $ getSubTreeUnsafe pth t)
       ]

-- Proposals for the rate tree.
proposalsRateTree :: Show a => Tree e a -> [Proposal I]
proposalsRateTree t =
  (rateTree @~ pulley 0.1 "Rate tree root" 5 True) :
  [ (rateTree . subTreeAt pth)
      @~ slideBranch 0.1 ("Rate tree branch " ++ show lb) 1 True
    | (pth, lb) <- itoList $ identify t,
      -- Since the stem does not change the likelihood, it is set to zero, and
      -- we do not slide the stem.
      not (null pth)
  ]
    ++ [ (rateTree . subTreeAt pth)
           @~ scaleTree 100 ("Rate tree node " ++ show lb) 1 True
         | (pth, lb) <- itoList $ identify t,
           -- Path does not lead to a leaf.
           not (null $ forest $ current $ goPathUnsafe pth $ fromTree t)
       ]

-- Create an accessor for a contrary proposal, see below.
timeHeightRateNormPair :: Lens' I (Double, Double)
timeHeightRateNormPair =
  lens
    (\x -> (x ^. timeHeight, x ^. rateMean))
    (\x (h, mu) -> x {_timeHeight = h, _rateMean = mu})

-- | The complete cycle includes proposals for the other parameters.
proposals :: Show a => Tree e a -> Cycle I
proposals t =
  fromList $
    [ timeBirthRate @~ scaleUnbiased 10 "Time birth rate" 10 True,
      timeDeathRate @~ scaleUnbiased 10 "Time death rate" 10 True,
      timeHeight @~ scaleUnbiased 3000 "Time height" 10 True,
      rateMean @~ scaleUnbiased 10 "Rate mean" 10 True,
      rateVariance @~ scaleUnbiased 10 "Rate variance" 10 True,
      timeHeightRateNormPair @~ scaleContrarily 10 0.1 "Time height, rate mean" 10 True
    ]
      ++ proposalsTimeTree t
      ++ proposalsRateTree t

-- Monitor the average rate. Useful, because it should not deviate from 1.0 too
-- much.
getRateAverage :: I -> Double
getRateAverage x = (/ 353) $ totalBranchLength $ x ^. rateTree

monParams :: [MonitorParameter I]
monParams =
  [ _timeBirthRate >$< monitorDouble "TimeBirthRate",
    _timeDeathRate >$< monitorDouble "TimeDeathRate",
    _timeHeight >$< monitorDouble "TimeHeight",
    _rateMean >$< monitorDouble "RateMean",
    _rateVariance >$< monitorDouble "RateVariance",
    getRateAverage >$< monitorDouble "RateAverage"
  ]

monStdOut :: MonitorStdOut I
monStdOut = monitorStdOut monParams 1

getTimeTreeNodeHeight :: Path -> I -> Double
getTimeTreeNodeHeight p x = (* h) $ getHeight $ label $ getSubTreeUnsafe p t
  where
    t = x ^. timeTree
    h = x ^. timeHeight

monCalibratedNodes :: [Calibration] -> [MonitorParameter I]
monCalibratedNodes cb = [getTimeTreeNodeHeight p >$< monitorDouble (name n a b) | (n, p, a, b) <- cb]
  where
    name s l r = "Calibration " ++ s ++ " (" ++ show l ++ ", " ++ show r ++ ")"

-- Positive if constraint is honored.
getTimeTreeDeltaNodeHeight :: Path -> Path -> I -> Double
getTimeTreeDeltaNodeHeight y o x = getTimeTreeNodeHeight o x - getTimeTreeNodeHeight y x

monConstrainedNodes :: [Constraint] -> [MonitorParameter I]
monConstrainedNodes cs = [getTimeTreeDeltaNodeHeight y o >$< monitorDouble (name n) | (n, y, o) <- cs]
  where
    name s = "Constraint " ++ s

monFileParams :: [Calibration] -> [Constraint] -> MonitorFile I
monFileParams cb cs =
  monitorFile
    "-params"
    ( monParams
        ++ monCalibratedNodes cb
        ++ monConstrainedNodes cs
    )
    1

absoluteTimeTree :: I -> Tree Double BS.ByteString
absoluteTimeTree s = first (* h) $ fromTimeTree t
  where
    h = s ^. timeHeight
    t = s ^. timeTree

monFileTimeTree :: MonitorFile I
monFileTimeTree = monitorFile "-timetree" [absoluteTimeTree >$< monitorTree "TimeTree"] 1

monFileRateTree :: MonitorFile I
monFileRateTree = monitorFile "-ratetree" [_rateTree >$< monitorTree "RateTree"] 1

-- | Monitor to standard output and files, as well as batch monitors.
monitor :: [Calibration] -> [Constraint] -> Monitor I
monitor cb cs = Monitor monStdOut [monFileParams cb cs, monFileTimeTree, monFileRateTree] []

-- | Number of burn in iterations.
nBurnIn :: Maybe Int
-- nBurnIn = Just 30
nBurnIn = Just 3000

-- | Auto tuning period.
nAutoTune :: Maybe Int
-- nAutoTune = Just 10
nAutoTune = Just 100

-- | Number of Metropolis-Hasting iterations after burn in.
nIterations :: Int
-- nIterations = 10
nIterations = 40000
