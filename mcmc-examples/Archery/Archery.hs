-- |
-- Module      :  Main
-- Description :  Modeling shots of an archer shots on a circular target
-- Copyright   :  (c) Dominik Schrempf 2020
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
import Numeric.Log
import qualified Statistics.Distribution as S
import qualified Statistics.Distribution.Exponential as S
import System.Random.MWC

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
distances :: GenIO -> IO [Distance]
distances = replicateM nArrows . S.genContVar (S.exponential muTrue)

-- Uninformative, improper prior for positive precision values.
pr :: PriorFunction Precision
pr = positive

-- Likelihood function.
lh :: [Distance] -> LikelihoodFunction Precision
lh xs p
  | p <= 0 = 0
  | otherwise = product [exponential p x | x <- xs]

-- The proposal cycle consists of one proposal only. A uniform distribution is used to
-- slide the precision of the archer.
proposals :: Cycle Precision
proposals = fromList [slideUniformSymmetric 1.0 (PName "Mu") (PWeight 1) Tune]

-- proposals = fromList [scaleUnbiased 1.6 (PName "Mu") (PWeight 1) Tune]
-- proposals = fromList [slide 0.06 0.8 (PName "Mu") (PWeight 1) Tune]
-- proposals = fromList [scaleBactrian 0.3 0.3 (PName "Mu") (PWeight 1) Tune]

-- Monitor the precision of the archer.
monMu :: MonitorParameter Precision
monMu = monitorDouble "Mu"

-- Monitor to standard output.
monStd :: MonitorStdOut Precision
monStd = monitorStdOut [monMu] 5000

-- Monitor to file.
monFile :: MonitorFile Precision
monFile = monitorFile "-mu" [monMu] 500

-- Monitor the batch mean of the precision of the archer.
monMuBatch :: MonitorParameterBatch Precision
monMuBatch = monitorBatchMean "Mean mu"

-- Monitor the batch mean to file.
monBatch :: MonitorBatch Precision
monBatch = monitorBatch "-mu" [monMuBatch] 1000

-- Combine the monitors.
mon :: Monitor Precision
mon = Monitor monStd [monFile] [monBatch]

-- Number of burn in iterations.
burnInSpec :: BurnIn
burnInSpec = BurnInWithAutoTuning 200000 10000

-- Number of Metropolis-Hastings iterations after burn in.
nIter :: Int
nIter = 1000000

main :: IO ()
main = do
  g <- create
  -- Simulate a list of observed arrow distances.
  xs <- distances g
  -- Combine all the objects defined above.
  let s = Settings "archery" burnInSpec nIter Overwrite NoSave Info
      c = chain pr (lh xs) proposals mon 0.01 g
  -- Run the Markov chain Monte Carlo sampler using the Metropolis-Hastings algorithm.
  void $ mcmcWith s (MHG c)
