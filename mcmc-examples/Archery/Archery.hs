-- |
-- Module      :  Main
-- Description :  Modeling shots of an archer shots on a circular target
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 27 17:11:28 2020.
--
-- Please refer to the corresponding [RevBayes
-- tutorial](https://revbayes.github.io/tutorials/mcmc/archery.html). Here, the
-- gamma distribution is not used to reduce computation of the posterior.
module Main
  ( main,
  )
where

import Control.Monad
import Mcmc
import qualified Statistics.Distribution as S
import qualified Statistics.Distribution.Exponential as S
import System.Random.Stateful

-- State space of the Markov chain. The precision of the archer is measured as a
-- 'Double'.
type Precision = Double

-- Distance of an arrow to the center.
type Distance = Double

-- Observed number of arrows.
nArrows :: Int
nArrows = 200

-- True precision of the archer.
muTrue :: Double
muTrue = 1.0

-- Simulated distances from center.
distances :: StatefulGen g m => g -> m [Distance]
distances = replicateM nArrows . S.genContVar (S.exponential muTrue)

-- Uninformative, improper prior for positive precision values.
pr :: PriorFunction Precision
pr = exponential 1.0

-- Likelihood function.
lh :: [Distance] -> LikelihoodFunction Precision
lh xs p
  | p <= 0 = 0
  | otherwise = product [exponential p x | x <- xs]

-- The proposal cycle consists of one proposal only. A uniform distribution is used to
-- slide the precision of the archer.
cc :: Cycle Precision
cc = cycleFromList [slideUniformSymmetric 1.0 (PName "Mu") (pWeight 1) Tune]

-- -- Other possibilities for proposals.
-- proposals = cycleFromList [scaleUnbiased 1.6 (PName "Mu") (PWeight 1) Tune]
-- proposals = cycleFromList [slide 0.06 0.8 (PName "Mu") (PWeight 1) Tune]
-- proposals = cycleFromList [scaleBactrian 0.3 0.3 (PName "Mu") (PWeight 1) Tune]

-- Monitor the precision of the archer.
monMu :: MonitorParameter Precision
monMu = monitorDouble "Mu"

-- Monitor to standard output.
monStd :: MonitorStdOut Precision
monStd = monitorStdOut [monMu] 5000

-- Monitor to file.
monFile :: MonitorFile Precision
monFile = monitorFile "mu" [monMu] 500

-- Monitor the batch mean of the precision of the archer.
monMuBatch :: MonitorParameterBatch Precision
monMuBatch = monitorBatchMean "Mean mu"

-- Monitor the batch mean to file.
monBatch :: MonitorBatch Precision
monBatch = monitorBatch "mu" [monMuBatch] 1000

-- Combine the monitors.
mon :: Monitor Precision
mon = Monitor monStd [monFile] [monBatch]

main :: IO ()
main = do
  let g = mkStdGen 0
  -- Simulate a list of observed arrow distances.
  let xs = runStateGen_ g distances
  -- MCMC settings and algorithm.
  let s =
        Settings
          (AnalysisName "archery")
          (BurnInWithAutoTuning 200000 10000)
          (Iterations 1000000)
          TraceAuto
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Info
  -- Use the Metropolis-Hastings-Green (MHG) algorithm.
  a <- mhg s pr (lh xs) cc mon 0.01 g
  -- Run the MCMC sampler.
  void $ mcmc s a
  -- Marginal likelihood estimation.
  let ss =
        MLSettings
          (AnalysisName "archery-marginal-likelihood")
          SteppingStoneSampling
          (NPoints 41)
          (BurnInWithAutoTuning 4000 100)
          (BurnInWithAutoTuning 1000 100)
          (Iterations 6000)
          Overwrite
          LogStdOutAndFile
          Info
  void $ marginalLikelihood ss pr (lh xs) cc (simpleMonitor 1000) 0.01 g
