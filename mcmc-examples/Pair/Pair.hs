-- |
-- Module      :  Main
-- Description :  Simple tests
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- This example module is more involved and contains proposals with custom
-- Jacobians. We sample the posterior distribution of a variable \(z \sim
-- \mbox{Exp}(1.0)\) separated into a sum \(z=x+y\), with \(x,y > 0\). The
-- sampler propose new values for \(x\) and \(y\).
module Main
  ( main,
  )
where

import Control.Lens
import Control.Monad
import Mcmc
import System.Random.MWC hiding (uniform)

-- The state is composed of a tuple of two variables x and y.
type I = (Double, Double)

analysisName :: AnalysisName
analysisName = AnalysisName "Pair"

-- Improper prior for positive values.
pr :: I -> Log Double
pr (x, y) = largerThan 0.00001 x * largerThan 0.00001 y

-- The likelihood function only acts on the sum of x and y, so we will need a
-- custom Jacobian in our proposals.
lh :: I -> Log Double
lh (x, y) = exponential 1.0 (x + y)

-- Initial value.
start :: I
start = (1.1, 1.1)

-- The Jacobian function required
jacobian :: I -> Log Double
jacobian = Exp . log . recip . uncurry (+)

-- The proposal cycle consists of one proposal only. A uniform distribution is used to
-- slide the precision of the archer.
cc :: Cycle I
cc =
  cycleFromList
    [
      -- The proposals require a custom Jacobian.
      liftProposalWith jacobian _1 (scaleUnbiased 1.0 (PName "x") (PWeight 1) Tune),
      liftProposalWith jacobian _2 (scaleUnbiased 1.0 (PName "y") (PWeight 1) Tune),
      -- Sliding the pair contrarily will not change z, and so, no Jacobian is
      -- required. If 'scaleContrarily' is used, a custom Jacobian is required.
      slideContrarily 0 0.5 (PName "x y") (PWeight 20) Tune
    ]

-- Monitor the precision of the archer.
monPs :: [MonitorParameter I]
monPs =
  [ fst >$< monitorDouble "x",
    snd >$< monitorDouble "y",
    uncurry (+) >$< monitorDouble "x+y",
    uncurry (*) >$< monitorDouble "x*y"
  ]

-- Monitor to standard output.
monStd :: MonitorStdOut I
monStd = monitorStdOut monPs 1000

monFile :: MonitorFile I
monFile = monitorFile "" monPs 100

-- Combine the monitors.
mon :: Monitor I
mon = Monitor monStd [monFile] []

main :: IO ()
main = do
  g <- create
  let mcmcS =
        Settings
          analysisName
          (BurnInWithAutoTuning 100000 500)
          (Iterations 500000)
          Overwrite
          Sequential
          NoSave
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg pr lh cc mon TraceAuto start g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings 3 1
  -- a <- mc3 mc3S pr lh cc mon start g
  void $ mcmc mcmcS a
