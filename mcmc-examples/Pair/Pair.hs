-- |
-- Module      :  Main
-- Description :  Simple tests
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
module Main
  ( main,
  )
where

import Control.Lens
import Control.Monad
import Mcmc
import System.Random.MWC hiding (uniform)

type I = (Double, Double)

-- Uninformative, improper prior for positive precision values.
pr :: I -> Log Double
pr (x, y) = largerThan 0.001 x * largerThan 0.001 y * exponential 1.0 (x * y)

-- Likelihood function.
lh :: I -> Log Double
lh _ = 1.0

-- The proposal cycle consists of one proposal only. A uniform distribution is used to
-- slide the precision of the archer.
cc :: Cycle I
cc =
  cycleFromList
    [ _1 @~ slideSymmetric 1 (PName "x") (PWeight 5) Tune,
      _2 @~ slideSymmetric 1 (PName "y") (PWeight 5) Tune,
      scaleContrarily 1.0 1.0 (PName "x y") (PWeight 5) Tune
    ]

-- Monitor the precision of the archer.
monPs :: [MonitorParameter I]
monPs =
  [ fst >$< monitorDouble "x",
    snd >$< monitorDouble "y",
    uncurry (*) >$< monitorDouble "xy"
  ]

-- Monitor to standard output.
monStd :: MonitorStdOut I
monStd = monitorStdOut monPs 1000

monFile :: MonitorFile I
monFile = monitorFile "" monPs 50

-- Combine the monitors.
mon :: Monitor I
mon = Monitor monStd [monFile] []

main :: IO ()
main = do
  g <- create
  let mcmcS =
        Settings
          (AnalysisName "test-pair")
          (BurnInWithAutoTuning 20000 1000)
          (Iterations 100000)
          Overwrite
          Sequential
          NoSave
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg pr lh cc mon TraceAuto (1, 1) g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings 3 1
  -- a <- mc3 mc3S pr lh cc mon (1, 1) g
  void $ mcmc mcmcS a
