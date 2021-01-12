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

-- TODO: Stepping stone sampling.
--
-- See Xie2010 and Fan2010.

import Control.Monad.IO.Class
import qualified Control.Monad.Parallel as P
import Control.Monad.Trans.Reader
import Data.Aeson
import qualified Data.Vector as VB
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
import Numeric.Log
import System.Random.MWC
import Text.Printf
import Prelude hiding (cycle)

-- Reciprocal temperature value traversed along the path integral.
type Point = Log Double

-- List of reciprocal temperature values traversed along the path integral.
type Points = [Point]

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

-- Distribute the points according to a skewed beta distribution (more points at
-- low values). It is inconvenient that the reciprocal temperatures are denoted
-- as beta, and we also use the beta distribution :). Don't mix them up!
--
-- See discussion in Xie, W., Lewis, P. O., Fan, Y., Kuo, L., & Chen, M.,
-- Improving marginal likelihood estimation for bayesian phylogenetic model
-- selection, Systematic Biology, 60(2), 150–160 (2010).
-- http://dx.doi.org/10.1093/sysbio/syq085
--
-- Or Figure 1 in Höhna, S., Landis, M. J., & Huelsenbeck, J. P., Parallel power
-- posterior analyses for fast computation of marginal likelihoods in
-- phylogenetics (2017). http://dx.doi.org/10.1101/104422
getPoints :: NPoints -> Points
getPoints x = [f i ** (1.0 / 0.3) | i <- [0 .. k1]]
  where
    k = fromNPoints x
    k1 = pred k
    f j = Exp $ log $ fromIntegral j / fromIntegral k1

-- TODO: Check acceptance ratio and warn if low or high.
sampleAtPoint ::
  ToJSON a =>
  Point ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  ML (MHG a)
sampleAtPoint b ss lhf a = do
  liftIO $ mcmc ss' a'
  where
    -- For debugging set a proper analysis name.
    nm = sAnalysisName ss
    getName :: Point -> AnalysisName
    getName x = nm <> AnalysisName (printf "%.5f" $ exp $ ln x)
    ss' = ss {sAnalysisName = getName b}
    -- Amend the likelihood function. Don't calculate the likelihood when beta
    -- is 0.0.
    lhf' = if b == 0.0 then const 1.0 else (** b) . lhf
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
  Points ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  -- Posterior probabilities.
  ML [Log Double]
traversePoints [] _ _ _ = return []
traversePoints (b : bs) ss lhf a = do
  -- Go to the next point.
  a' <- sampleAtPoint b ss lhf a
  -- Extract the links.
  ls <- liftIO $ takeT n $ trace $ fromMHG a'
  -- Calculate the mean posterior probability at the point.
  let mp = getMeanPosterior ls
  -- Get the mean posterior probabilities of the other points.
  mps <- traversePoints bs ss lhf a'
  return $ mp : mps
  where
    n = fromIterations $ sIterations ss
    getPosterior x = prior x * likelihood x
    getMeanPosterior xs = VB.sum (VB.map getPosterior xs) / fromIntegral (VB.length xs)

-- TODO: Proper return value; marginal likelihood and confidence interval.

-- TODO: Proper output.

mlTIRun ::
  ToJSON a =>
  Points ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  a ->
  GenIO ->
  -- Mean posteriors likelihoods at points.
  ML [Log Double]
mlTIRun xs prf lhf cc i0 g = do
  s <- reader settings
  let nm = mlAnalysisName s
      is = mlIterations s
      biI = mlInitialBurnIn s
      -- TODO: Maybe take the execution mode and verbosity from MLSettings.
      ssI = Settings nm biI is Fail Sequential NoSave Quiet
      biP = mlPointBurnIn s
      -- TODO: Maybe take the execution mode and verbosity from MLSettings.
      ssP = Settings nm biP is Fail Sequential NoSave Quiet
      trLen = TraceMinimum $ fromIterations is
      mn = noMonitor 1
  a0 <- liftIO $ mhg prf lhf cc mn trLen i0 g
  a1 <- sampleAtPoint x0 ssI lhf a0
  traversePoints xs ssP lhf a1
  where
    x0 = head xs

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
  -- | Marginal likelihoods of the path integrals in forward and backward
  -- direction.
  IO (Log Double, Log Double)
mlThermodynamicIntegration s prf lhf cc i0 g = do
  -- Initialize.
  e <- initializeEnvironment s
  [g0, g1] <- splitGen 2 g

  -- Run.
  runReaderT
    ( do
        -- Parallel execution of both path integrals.
        [mpsForward, mpsBackward] <-
          P.sequence
            [ mlTIRun bsForward prf lhf cc i0 g0,
              mlTIRun bsBackward prf lhf cc i0 g1
            ]
        -- We have to calculate the integral here, 'triangle' chokes when the
        -- points are reversed (log domain values cannot handle negative
        -- values).
        return
          ( triangle mpsForward bsForward,
            triangle (reverse mpsBackward) bsForward
          )
    )
    e
  where
    bsForward = getPoints $ mlNPoints s
    bsBackward = reverse bsForward

triangle ::
  -- Y values.
  [Log Double] ->
  -- X values.
  Points ->
  -- Integral.
  Log Double
triangle (x0 : x1 : xs) (b0 : b1 : bs) = (x0 + x1) / (b1 - b0) + triangle (x1 : xs) (b1 : bs)
triangle _ _ = 0

-- TODO.

-- | Stub.
mlSteppingStoneSampling :: IO ()
mlSteppingStoneSampling = undefined
