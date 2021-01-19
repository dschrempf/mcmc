{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.MarginalLikelihood
-- Description :  Calculate the marginal likelihood
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jan 11 16:34:18 2021.
module Mcmc.MarginalLikelihood
  ( NPoints (..),
    MLSettings (..),
    mlThermodynamicIntegration,
    mlSteppingStoneSampling,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import qualified Control.Monad.Parallel as P
import Control.Monad.Trans.Reader
import Data.Aeson
import qualified Data.Vector as VB
import qualified Data.Vector.Unboxed as VU
import Mcmc.Algorithm.MHG
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Trace
import Mcmc.Environment
import Mcmc.Internal.Random
import Mcmc.Logger
import Mcmc.Mcmc
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log hiding (sum)
import System.Random.MWC
import Text.Printf
import Text.Show.Pretty
import Prelude hiding (cycle)

-- Reciprocal temperature value traversed along the path integral.
type Point = Double

-- | The number of points used to approximate the path integral.
newtype NPoints = NPoints {fromNPoints :: Int}
  deriving (Eq, Read, Show)

data MLSettings = MLSettings
  { mlAnalysisName :: AnalysisName,
    mlNPoints :: NPoints,
    -- | Initial burn in at the starting point of the path.
    mlInitialBurnIn :: BurnInSpecification,
    -- | Repetitive burn in at each point on the path.
    mlPointBurnIn :: BurnInSpecification,
    -- | The number of iterations performed at each point.
    mlIterations :: Iterations,
    mlExecutionMode :: ExecutionMode,
    mlVerbosity :: Verbosity
  }
  deriving (Eq, Read, Show)

instance HasAnalysisName MLSettings where
  getAnalysisName = mlAnalysisName

instance HasExecutionMode MLSettings where
  getExecutionMode = mlExecutionMode

instance HasVerbosity MLSettings where
  getVerbosity = mlVerbosity

type ML a = ReaderT (Environment MLSettings) IO a

-- See 'getPoints'. Alpha=0.3 is the standard choice.
alpha :: Double
alpha = 0.3

-- Distribute the points according to a skewed beta distribution with given
-- 'alpha' value. If alpha is below 1.0, more points at lower values, which is
-- desired. It is inconvenient that the reciprocal temperatures are denoted as
-- beta, and we also use the beta distribution :). Don't mix them up!
--
-- See discussion in Xie, W., Lewis, P. O., Fan, Y., Kuo, L., & Chen, M.,
-- Improving marginal likelihood estimation for bayesian phylogenetic model
-- selection, Systematic Biology, 60(2), 150–160 (2010).
-- http://dx.doi.org/10.1093/sysbio/syq085
--
-- Or Figure 1 in Höhna, S., Landis, M. J., & Huelsenbeck, J. P., Parallel power
-- posterior analyses for fast computation of marginal likelihoods in
-- phylogenetics (2017). http://dx.doi.org/10.1101/104422
getPoints :: NPoints -> [Point]
getPoints x = [f i ** (1.0 / alpha) | i <- [0 .. k1]]
  where
    k = fromNPoints x
    k1 = pred k
    f j = fromIntegral j / fromIntegral k1

sampleAtPoint ::
  ToJSON a =>
  Point ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  ML (MHG a)
sampleAtPoint x ss lhf a = do
  a'' <- liftIO $ mcmc ss' a'
  let ar = acceptanceRates $ acceptance $ fromMHG a''
      meanAr = sum ar / fromIntegral (length ar)
  when (meanAr <= 0.1) $ logWarnB "Overall acceptance rate is below 0.1."
  when (meanAr >= 0.9) $ logWarnB "Overall acceptance rate is above 0.9."
  return a''
  where
    -- For debugging set a proper analysis name.
    nm = sAnalysisName ss
    getName :: Point -> AnalysisName
    getName y = nm <> AnalysisName (printf "%.5f" y)
    ss' = ss {sAnalysisName = getName x}
    -- Amend the likelihood function. Don't calculate the likelihood when the
    -- point is 0.0.
    lhf' = if x == 0.0 then const 1.0 else (** Exp (log x)) . lhf
    -- Amend the MHG algorithm.
    ch = fromMHG a
    l = link ch
    ch' =
      ch
        { -- Important: Update the likelihood using the new likelihood function.
          link = l {likelihood = lhf' $ state l},
          iteration = 0,
          start = 0,
          likelihoodFunction = lhf'
        }
    a' = MHG ch'

traversePoints ::
  ToJSON a =>
  [Point] ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  -- For each point a vector of obtained likelihoods.
  ML [VU.Vector (Log Double)]
traversePoints [] _ _ _ = return []
traversePoints (b : bs) ss lhf a = do
  -- Sample the next point.
  a' <- sampleAtPoint b ss lhf a
  -- Get the links samples at this point.
  ls <- liftIO $ takeT n $ trace $ fromMHG a'
  -- Extract the likelihoods.
  let lhs = VU.convert $ VB.map (lhf . state) ls
  -- Sample the other points.
  lhss <- traversePoints bs ss lhf a'
  return $ lhs : lhss
  where
    n = fromIterations $ sIterations ss

mlTIRun ::
  ToJSON a =>
  [Point] ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  a ->
  GenIO ->
  -- Mean log likelihoods at points.
  ML [Double]
mlTIRun xs prf lhf cc i0 g = do
  logDebugB "mlTiRun: Begin."
  s <- reader settings
  let nm = mlAnalysisName s
      is = mlIterations s
      biI = mlInitialBurnIn s
      biP = mlPointBurnIn s
      -- Be quiet for the sub MCMC runs.
      ssI = Settings nm biI is Fail Sequential NoSave Quiet
      ssP = Settings nm biP is Fail Sequential NoSave Quiet
      trLen = TraceMinimum $ fromIterations is
      mn = noMonitor 1
  logDebugB "mlTiRun: Initialize MHG algorithm."
  a0 <- liftIO $ mhg prf lhf cc mn trLen i0 g
  logDebugB "mlTiRun: Sample first point."
  a1 <- sampleAtPoint x0 ssI lhf a0
  logDebugB "mlTiRun: Traverse points."
  lhss <- traversePoints xs ssP lhf a1
  logDebugB "mlTiRun: Calculate mean log likelihoods."
  -- It is important to average across the log likelihoods here (and not the
  -- likelihoods). I am not exactly sure why this is.
  let meanLogLhs = map (\x -> VU.sum (VU.map ln x) / fromIntegral (VU.length x)) lhss
  return meanLogLhs
  where
    x0 = head xs

-- TODO: Confidence interval.

-- | Calculate the marginal likelihood using a path integral.
--
-- Also known as thermodynamic integration. In particular, /Annealing-Melting
-- Integration/ is used.
--
-- See Lartillot, N., & Philippe, H., Computing Bayes Factors Using
-- Thermodynamic Integration, Systematic Biology, 55(2), 195–207 (2006).
-- http://dx.doi.org/10.1080/10635150500433722
mlThermodynamicIntegration ::
  ToJSON a =>
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  -- | Initial state.
  a ->
  GenIO ->
  -- | Marginal log likelihoods of the path integrals in forward and backward
  -- direction.
  IO (Double, Double)
mlThermodynamicIntegration s prf lhf cc i0 g = do
  -- Initialize.
  e <- initializeEnvironment s
  [g0, g1] <- splitGen 2 g

  -- Run.
  runReaderT
    ( do
        logInfoStartingTime
        logInfoB "Estimate marginal likelihood using thermodynamic integration."
        logDebugB "The marginal likelihood settings are:"
        logDebugS $ ppShow s
        -- Parallel execution of both path integrals.
        [mllhsForward, mllhsBackward] <-
          P.sequence
            [ mlTIRun bsForward prf lhf cc i0 g0,
              mlTIRun bsBackward prf lhf cc i0 g1
            ]
        logInfoEndTime
        let mlForward = integrateSimpsonTriangle bsForward mllhsForward
            mlBackward = negate $ integrateSimpsonTriangle bsBackward mllhsBackward
        logInfoB "Marginal log likelihood:"
        logInfoS $ "Forward:  " ++ show mlForward
        logInfoS $ "Backward: " ++ show mlBackward
        return (mlForward, mlBackward)
    )
    e
  where
    bsForward = getPoints $ mlNPoints s
    bsBackward = reverse bsForward

-- Use lists since the number of points is expected to be low.
integrateSimpsonTriangle ::
  -- X values.
  [Point] ->
  -- Y values.
  [Double] ->
  -- Integral.
  Double
integrateSimpsonTriangle xs ys = 0.5 * go xs ys
  where
    go (p0 : p1 : ps) (z0 : z1 : zs) = (z0 + z1) * (p1 - p0) + go (p1 : ps) (z1 : zs)
    go _ _ = 0

-- TODO: See Xie2010 and Fan2010.

-- | Stub.
mlSteppingStoneSampling :: IO ()
mlSteppingStoneSampling = undefined
