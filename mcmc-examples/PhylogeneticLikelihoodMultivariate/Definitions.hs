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

import Control.Comonad
import Control.Lens
import Data.Aeson
import Data.Bifunctor
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
-- The topologies of the time and rate tree are equal. This is, however, not
-- ensured by the types. In the future, we may just use one tree storing both,
-- the times and the rates.
--
-- Remark:
--
-- Let T be the length of a branch measured in unit time, and R be the rate on
-- this branch. Then, the length measured in number of substitutions of this
-- very same branch is d=T*R. Internally, a relative time t and relative rate r
-- are used to measure the branch length such that d=T*R=(T/h)*(h*R):=t*r, where
-- h is the root height of the tree measured in unit time.
--
-- In brief, the relative time and rate are defined as t=T/h, and r=R*h, where h
-- is the root height.
--
-- This has two big advantages:
--
-- 1. The time tree object storing the relative times is a normalized tree with
--    root height 1.0.
--
-- 2. The likelihood can be calculated without consulting the root height h
--    measured in unit time. This is important, because there is simply no
--    information about the root age in the alignment which is used to calculate
--    the likelihood.
--
-- In turn, the tree height h measured in unit time is only determined by the
-- calibrations, the constraints, and the birth and death process prior.
--
-- I think this is a clean solution.
data I = I
  { -- | Birth rate of time tree.
    _timeBirthRate :: Double,
    -- | Death rate of time tree.
    _timeDeathRate :: Double,
    -- | Height of the time tree in absolute time measured in million years, see
    -- calibrations.
    _timeHeight :: Double,
    -- | Normalized time tree of height 1.0. Branch labels denote relative
    -- times; node labels denote relative node height.
    _timeTree :: Tree Double Double,
    -- | Shape of the gamma distribution prior of the rates.
    _rateShape :: Double,
    -- | Scale of the gamma distribution prior of the rates.
    _rateScale :: Double,
    -- | Rate tree. Branch labels denote relative rates; node labels are unused.
    _rateTree :: Tree Double ()
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
cleanTimeTree = timeTree %~ (extend rootHeight . normalizeHeight . makeUltrametric)

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
initWith :: Tree Double Int -> I
initWith t =
  I
    { _timeBirthRate = 1.0,
      _timeDeathRate = 1.0,
      _timeHeight = 1200.0,
      _timeTree = t',
      _rateShape = 10.0,
      _rateScale = 2.0,
      _rateTree = bimap (const 1.0) (const ()) t
    }
  where
    t' = extend rootHeight $ normalizeHeight $ makeUltrametric t

-- Calibrations are defined in the module 'Calibration'.

-- Constraints are defined in the module 'Constraint'.

-- | Prior distribution.
priorDistribution :: [Calibration] -> [Constraint] -> I -> Log Double
priorDistribution cb cs (I l m h t k th r) =
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

-- | The complete cycle includes proposals for the other parameters.
proposals :: Show a => Tree e a -> Cycle I
proposals t =
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

-- TODO: Remove; debug.
monFileTimeTreeWithHeight :: MonitorFile I
monFileTimeTreeWithHeight = monitorFile "-timetree-height" [_timeTree >$< monitorTree "TimeTreeWithHeight"] 1

monFileRateTree :: MonitorFile I
monFileRateTree = monitorFile "-ratetree" [_rateTree >$< monitorTree "RateTree"] 1

-- | Monitor to standard output and files, as well as batch monitors.
monitor :: Monitor I
monitor = Monitor monStdOut [monFileParams, monFileTimeTree, monFileRateTree, monFileTimeTreeWithHeight] []

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
nIterations = 10000
