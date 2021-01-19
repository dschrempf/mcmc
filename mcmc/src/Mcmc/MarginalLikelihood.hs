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
    MLAlgorithm (..),
    MLSettings (..),
    marginalLikelihood,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import qualified Control.Monad.Parallel as P
import Control.Monad.Trans.Reader
import Data.Aeson
import Data.List
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

-- | Algorithms to calculate the marginal likelihood.
data MLAlgorithm
  = -- | Use a classical path integral. Also known as thermodynamic integration.
    -- In particular, /Annealing-Melting Integration/ is used.
    --
    -- See Lartillot, N., & Philippe, H., Computing Bayes Factors Using
    -- Thermodynamic Integration, Systematic Biology, 55(2), 195–207 (2006).
    -- http://dx.doi.org/10.1080/10635150500433722
    ThermodynamicIntegration
  | -- | Use stepping stone sampling.
    --
    -- See Xie, W., Lewis, P. O., Fan, Y., Kuo, L., & Chen, M., Improving
    -- marginal likelihood estimation for Bayesian phylogenetic model selection,
    -- Systematic Biology, 60(2), 150–160 (2010).
    -- http://dx.doi.org/10.1093/sysbio/syq085
    --
    -- Or Fan, Y., Wu, R., Chen, M., Kuo, L., & Lewis, P. O., Choosing among
    -- partition models in bayesian phylogenetics, Molecular Biology and
    -- Evolution, 28(1), 523–532 (2010). http://dx.doi.org/10.1093/molbev/msq224
    SteppingStoneSampling
  deriving (Eq, Read, Show)

data MLSettings = MLSettings
  { mlAnalysisName :: AnalysisName,
    mlAlgorithm :: MLAlgorithm,
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
  -- For each point a vector of obtained likelihoods stored in the log domain.
  ML [VU.Vector (Log Double)]
traversePoints [] _ _ _ = return []
traversePoints (b : bs) ss lhf a = do
  -- Sample the next point.
  a' <- sampleAtPoint b ss lhf a
  -- Get the links samples at this point.
  ls <- liftIO $ takeT n $ trace $ fromMHG a'
  -- Extract the likelihoods.
  --
  -- NOTE: This could be sped up by mapping (** -b) on the power likelihoods.
  let lhs = VU.convert $ VB.map (lhf . state) ls
  -- Sample the other points.
  lhss <- traversePoints bs ss lhf a'
  return $ lhs : lhss
  where
    n = fromIterations $ sIterations ss

mlRun ::
  ToJSON a =>
  [Point] ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  a ->
  GenIO ->
  -- For each point a vector of likelihoods stored in log domain.
  ML [VU.Vector (Log Double)]
mlRun xs prf lhf cc i0 g = do
  logDebugB "mlRun: Begin."
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
  logDebugB "mlRun: Initialize MHG algorithm."
  a0 <- liftIO $ mhg prf lhf cc mn trLen i0 g
  logDebugB "mlRun: Sample first point with initial burn in settings."
  a1 <- sampleAtPoint x0 ssI lhf a0
  logDebugB "mlRun: Traverse points."
  traversePoints xs ssP lhf a1
  where
    x0 = head xs

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

tiWrapper ::
  ToJSON a =>
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  a ->
  GenIO ->
  -- Marginal likelihood in log domain.
  ML (Log Double)
tiWrapper s prf lhf cc i0 g = do
  logInfoB "Path integral (thermodynamic integration)."
  [g0, g1] <- splitGen 2 g

  -- Parallel execution of both path integrals.
  [lhssForward, lhssBackward] <-
    P.sequence
      [ mlRun bsForward prf lhf cc i0 g0,
        mlRun bsBackward prf lhf cc i0 g1
      ]
  logInfoEndTime

  logDebugB "tiWrapper: Calculate mean log likelihoods."
  -- It is important to average across the log likelihoods here (and not the
  -- likelihoods). I am not exactly sure why this is.
  let getMeanLogLhs = map (\x -> VU.sum (VU.map ln x) / fromIntegral (VU.length x))
      mlForward = integrateSimpsonTriangle bsForward (getMeanLogLhs lhssForward)
      mlBackward = negate $ integrateSimpsonTriangle bsBackward (getMeanLogLhs lhssBackward)
  logDebugS $ "tiWrapper: Marginal log likelihood of forward integral: " ++ show mlForward
  logDebugS $ "tiWrapper: Marginal log likelihood of backward integral: " ++ show mlBackward
  let mean = 0.5 * (mlForward + mlBackward)
  logDebugS $ "tiWrapper: The mean is: " ++ show mean
  return $ Exp mean
  where
    bsForward = getPoints $ mlNPoints s
    bsBackward = reverse bsForward

-- Helper function to exponentiate log domain values with a double value.
pow' :: Log Double -> Double -> Log Double
pow' x p = Exp $ ln x * p

-- See Xie2010 p. 153, bottom left.
sssCalculateMarginalLikelihood :: [Point] -> [VU.Vector (Log Double)] -> Log Double
sssCalculateMarginalLikelihood xs lhss = product $ zipWith3 f xs (tail xs) lhss
  where
    f :: Point -> Point -> VU.Vector (Log Double) -> Log Double
    -- f beta_{k-1} beta_k lhs_{k-1}
    f bkm1 bk lhs = n1 * VU.sum lhsPowered
      where
        n1 = recip $ fromIntegral $ VU.length lhs
        dbeta = bk - bkm1
        lhsPowered = VU.map (`pow'` dbeta) lhs
    -- -- Numerical stability by factoring out lhMax. But no observed
    -- -- improvement towards the standard version.
    --
    -- f bkm1 bk lhs = n1 * pow' lhMax dbeta * VU.sum lhsNormedPowered
    --   where n1 = recip $ fromIntegral $ VU.length lhs
    --         lhMax = VU.maximum lhs
    --         dbeta = bk - bkm1
    --         lhsNormed = VU.map (/lhMax) lhs
    --         lhsNormedPowered = VU.map (`pow'` dbeta) lhsNormed

-- -- Computation of the log of the marginal likelihood. According to the paper,
-- -- this estimator is biased and I did not observe any improvements compared
-- -- to the direct estimator implemented above.
--
-- -- See Xie2010 p. 153, top right.
-- sssCalculateMarginalLikelihood' :: [Point] -> [VU.Vector (Log Double)] -> Log Double
-- sssCalculateMarginalLikelihood' xs lhss = Exp $ sum $ zipWith3 f xs (tail xs) lhss
--   where f :: Point -> Point -> VU.Vector (Log Double) -> Double
--         -- f beta_{k-1} beta_k lhs_{k-1}
--         f bkm1 bk lhs = dbeta * llhMax + log (n1 * VU.sum lhsNormedPowered)
--           where dbeta = bk - bkm1
--                 llhMax = ln $ VU.maximum lhs
--                 n1 = recip $ fromIntegral $ VU.length lhs
--                 llhs = VU.map ln lhs
--                 llhsNormed = VU.map (\x -> x - llhMax) llhs
--                 lhsNormedPowered = VU.map (\x -> exp $ dbeta * x) llhsNormed
sssWrapper ::
  ToJSON a =>
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  a ->
  GenIO ->
  -- Marginal likelihood in log domain.
  ML (Log Double)
sssWrapper s prf lhf cc i0 g = do
  logInfoB "Stepping stone sampling."
  -- The last point does not need to be sampled.
  logLhss <- mlRun bsForward' prf lhf cc i0 g
  logDebugB "sssWrapper: Calculate marginal likelihood."
  return $ sssCalculateMarginalLikelihood bsForward logLhss
  where
    bsForward = getPoints $ mlNPoints s
    bsForward' = init bsForward

marginalLikelihood ::
  ToJSON a =>
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  -- | Initial state.
  a ->
  GenIO ->
  -- | Estimated marginal likelihood stored in log domain.
  IO (Log Double)
marginalLikelihood s prf lhf cc i0 g = do
  -- Initialize.
  e <- initializeEnvironment s

  -- Run.
  runReaderT
    ( do
        logInfoStartingTime
        logInfoB "Estimate marginal likelihood."
        logDebugB "The marginal likelihood settings are:"
        logDebugS $ ppShow s
        val <- case mlAlgorithm s of
          ThermodynamicIntegration -> tiWrapper s prf lhf cc i0 g
          SteppingStoneSampling -> sssWrapper s prf lhf cc i0 g
        logInfoS $ "Marginal log likelihood: " ++ show (ln val)
        -- TODO: Simulation variance.
        logInfoS "The simulation variance is not yet available."
        return val
    )
    e
