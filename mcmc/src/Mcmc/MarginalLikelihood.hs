{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.MarginalLikelihood
-- Description :  Calculate the marginal likelihood
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jan 11 16:34:18 2021.
module Mcmc.MarginalLikelihood
  ( MarginalLikelihood,
    NPoints (..),
    MLAlgorithm (..),
    MLSettings (..),
    marginalLikelihood,
  )
where

import Control.Concurrent (getNumCapabilities)
import Control.Concurrent.Async hiding (link)
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
import Control.Monad.Trans.Reader
import Data.Aeson
import Data.List hiding (cycle)
import qualified Data.Map.Strict as M
import qualified Data.Vector as VB
import qualified Data.Vector.Unboxed as VU
import Mcmc.Acceptance
import Mcmc.Algorithm.MHG
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Trace
import Mcmc.Cycle
import Mcmc.Environment
import Mcmc.Likelihood
import Mcmc.Logger
import Mcmc.Mcmc
import Mcmc.Monitor
import Mcmc.Prior
import Mcmc.Settings
import Numeric.Log hiding (sum)
import System.Directory
import System.Random.Stateful
import Text.Printf
import Prelude hiding (cycle)

-- | Marginal likelihood values are stored in log domain.
type MarginalLikelihood = Log Double

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

-- | Settings of the marginal likelihood estimation.
data MLSettings = MLSettings
  { mlAnalysisName :: AnalysisName,
    mlAlgorithm :: MLAlgorithm,
    mlNPoints :: NPoints,
    -- | Initial burn in at the starting point of the path (or each segment if
    -- running in parallel).
    mlInitialBurnIn :: BurnInSettings,
    -- | Repetitive burn in at each point on the path.
    mlPointBurnIn :: BurnInSettings,
    -- | The number of iterations performed at each point.
    mlIterations :: Iterations,
    mlExecutionMode :: ExecutionMode,
    mlParallelizationMode :: ParallelizationMode,
    mlLogMode :: LogMode,
    mlVerbosity :: Verbosity
  }
  deriving (Eq, Read, Show)

instance HasAnalysisName MLSettings where
  getAnalysisName = mlAnalysisName

instance HasExecutionMode MLSettings where
  getExecutionMode = mlExecutionMode

instance HasLogMode MLSettings where
  getLogMode = mlLogMode

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
  (ToJSON a) =>
  Bool ->
  Point ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  ML (MHG a)
sampleAtPoint isInitialBurnIn x ss lhf a = do
  a'' <- liftIO $ mcmc ss' a'
  let ch'' = fromMHG a''
      ac = acceptances ch''
      mAr = sequence $ acceptanceRates ac
  logDebugB "sampleAtPoint: Summarize cycle."
  logDebugB $ summarizeCycle AllProposals ac $ cycle ch''
  unless
    isInitialBurnIn
    ( case mAr of
        Nothing -> logWarnB "Some acceptance rates are unavailable."
        Just ar -> do
          unless (M.null $ M.filter (<= 0.1) ar) $ logWarnB "Some acceptance rates are below 0.1."
          unless (M.null $ M.filter (>= 0.9) ar) $ logWarnB "Some acceptance rates are above 0.9."
    )
  return a''
  where
    -- For debugging set a proper analysis name.
    nm = sAnalysisName ss
    getName :: Point -> AnalysisName
    getName y = nm <> AnalysisName ("/" <> printf "point%.8f" y)
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
  (ToJSON a) =>
  NPoints ->
  [(Int, Point)] ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  -- For each point a vector of obtained likelihoods stored in the log domain.
  ML [VU.Vector Likelihood]
traversePoints _ [] _ _ _ = return []
traversePoints k ((idb, b) : bs) ss lhf a = do
  let msg = printf "Point %4d of %4d: %.12f." idb k' b
  logInfoS msg
  a' <- sampleAtPoint False b ss lhf a
  -- Get the links samples at this point.
  ls <- liftIO $ takeT n $ trace $ fromMHG a'
  -- Extract the likelihoods.
  --
  -- NOTE: This could be sped up by mapping (** -b) on the power likelihoods.
  --
  -- NOTE: This bang is an important one, because if the lhs are not strictly
  -- calculated here, the complete MCMC runs are dragged along before doing so
  -- resulting in a severe memory leak.
  let !lhs = VU.convert $ VB.map (lhf . state) ls
  -- Sample the other points.
  lhss <- traversePoints k bs ss lhf a'
  return $ lhs : lhss
  where
    n = fromIterations $ sIterations ss
    (NPoints k') = k

nChunks :: Int -> [a] -> [[a]]
nChunks k xs = chop (chunks k l) xs
  where
    l = length xs

chunks :: Int -> Int -> [Int]
chunks c n = filter (> 0) ns
  where
    n' = n `div` c
    r = n `mod` c
    ns = replicate r (n' + 1) ++ replicate (c - r) n'

chop :: [Int] -> [a] -> [[a]]
chop [] [] = []
chop (n : ns) xs
  | n > 0 = take n xs : chop ns (drop n xs)
  | otherwise = error "chop: n negative or zero"
chop _ _ = error "chop: not all list elements handled"

mlRunPar ::
  (ToJSON a) =>
  ParallelizationMode ->
  NPoints ->
  [(Int, Point)] ->
  ExecutionMode ->
  Verbosity ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  a ->
  StdGen ->
  ML [VU.Vector Likelihood]
mlRunPar pm k xs em vb prf lhf cc mn i0 g = do
  nThreads <- case pm of
    Sequential -> do
      logInfoB "Sequential execution."
      pure 1
    Parallel -> do
      n <- liftIO getNumCapabilities
      logInfoS $ "Parallel execution with " <> show n <> " cores."
      pure n
  let xsChunks = nChunks nThreads xs
  r <- ask
  xss <-
    liftIO $
      mapConcurrently
        (\thesePoints -> runReaderT (mlRun k thesePoints em vb prf lhf cc mn i0 g) r)
        xsChunks
  pure $ concat xss

mlRun ::
  (ToJSON a) =>
  NPoints ->
  [(Int, Point)] ->
  ExecutionMode ->
  Verbosity ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  a ->
  StdGen ->
  -- For each point a vector of likelihoods stored in log domain.
  ML [VU.Vector Likelihood]
mlRun k xs em vb prf lhf cc mn i0 g = do
  logDebugB "mlRun: Begin."
  s <- reader settings
  let nm = mlAnalysisName s
      is = mlIterations s
      biI = mlInitialBurnIn s
      biP = mlPointBurnIn s
      -- Only log sub MCMC samplers when debugging.
      vb' = case vb of
        Debug -> Debug
        _ -> Quiet
      trLen = TraceMinimum $ fromIterations is
      ssI = Settings nm biI (Iterations 0) trLen em Sequential NoSave LogFileOnly vb'
      ssP = Settings nm biP is trLen em Sequential NoSave LogFileOnly vb'
  logDebugB "mlRun: Initialize MHG algorithm."
  a0 <- liftIO $ mhg ssI prf lhf cc mn i0 g
  let msg = printf "Initial burn in at point %.12f with ID %4d." x0 id0
  logInfoS msg
  a1 <- sampleAtPoint True x0 ssI lhf a0
  logDebugB "mlRun: Traverse points."
  traversePoints k xs ssP lhf a1
  where
    (id0, x0) = head xs

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
  (ToJSON a) =>
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  a ->
  StdGen ->
  ML MarginalLikelihood
tiWrapper s prf lhf cc mn i0 g = do
  logInfoB "Path integral (thermodynamic integration)."
  let (g0, g1) = split g

  -- Parallel execution of both path integrals.
  r <- ask
  (lhssForward, lhssBackward) <-
    lift $
      concurrently
        (runReaderT (mlRunPar pm k (zip [1 ..] bsForward) em vb prf lhf cc mn i0 g0) r)
        (runReaderT (mlRunPar pm k (zip [1 ..] bsBackward) em vb prf lhf cc mn i0 g1) r)
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
    k = mlNPoints s
    bsForward = getPoints k
    bsBackward = reverse bsForward
    em = mlExecutionMode s
    pm = mlParallelizationMode s
    vb = mlVerbosity s

-- Helper function to exponentiate log domain values with a double value.
pow' :: Log Double -> Double -> Log Double
pow' x p = Exp $ ln x * p

-- See Xie2010 p. 153, bottom left.
sssCalculateMarginalLikelihood :: [Point] -> [VU.Vector Likelihood] -> MarginalLikelihood
sssCalculateMarginalLikelihood xs lhss = product $ zipWith3 f xs (tail xs) lhss
  where
    f :: Point -> Point -> VU.Vector Likelihood -> MarginalLikelihood
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
-- sssCalculateMarginalLikelihood' :: [Point] -> [VU.Vector Likelihood] -> MarginalLikelihood
-- sssCalculateMarginalLikelihood' xs lhss = Exp $ sum $ zipWith3 f xs (tail xs) lhss
--   where f :: Point -> Point -> VU.Vector Likelihood -> Double
--         -- f beta_{k-1} beta_k lhs_{k-1}
--         f bkm1 bk lhs = dbeta * llhMax + log (n1 * VU.sum lhsNormedPowered)
--           where dbeta = bk - bkm1
--                 llhMax = ln $ VU.maximum lhs
--                 n1 = recip $ fromIntegral $ VU.length lhs
--                 llhs = VU.map ln lhs
--                 llhsNormed = VU.map (\x -> x - llhMax) llhs
--                 lhsNormedPowered = VU.map (\x -> exp $ dbeta * x) llhsNormed
sssWrapper ::
  (ToJSON a) =>
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  a ->
  StdGen ->
  ML MarginalLikelihood
sssWrapper s prf lhf cc mn i0 g = do
  logInfoB "Stepping stone sampling."
  logLhss <- mlRunPar pm k (zip [1 ..] bsForward') em vb prf lhf cc mn i0 g
  logInfoB "The last point does not need to be sampled with stepping stone sampling."
  logDebugB "sssWrapper: Calculate marginal likelihood."
  return $ sssCalculateMarginalLikelihood bsForward logLhss
  where
    k = mlNPoints s
    bsForward = getPoints k
    bsForward' = init bsForward
    em = mlExecutionMode s
    pm = mlParallelizationMode s
    vb = mlVerbosity s

-- | Estimate the marginal likelihood.
marginalLikelihood ::
  (ToJSON a) =>
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  InitialState a ->
  StdGen ->
  IO MarginalLikelihood
marginalLikelihood s prf lhf cc mn i0 g = do
  -- Initialize.
  e <- initializeEnvironment s

  -- Create marginal likelihood analysis directory.
  let n = fromAnalysisName $ mlAnalysisName s
  createDirectoryIfMissing True n

  -- Run.
  runReaderT
    ( do
        logInfoStartingTime
        logInfoB "Estimate marginal likelihood."
        logDebugB "marginalLikelihood: The marginal likelihood settings are:"
        logDebugS $ show s
        val <- case mlAlgorithm s of
          ThermodynamicIntegration -> tiWrapper s prf lhf cc mn i0 g
          SteppingStoneSampling -> sssWrapper s prf lhf cc mn i0 g
        logInfoS $ "Marginal log likelihood: " ++ show (ln val)
        -- TODO (low): Simulation variance.
        logInfoS "The simulation variance is not yet available."
        return val
    )
    e
