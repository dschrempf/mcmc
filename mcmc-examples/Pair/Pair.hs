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
-- This example module is more involved and contains custom proposals with
-- custom Jacobians. We sample the posterior distribution of a variable \(z \sim
-- \mbox{Exp}(1.0)\) separated into a sum \(z=x+y\), with \(x,y > 0\). The
-- sampler propose new values for \(x\) and \(y\).
module Main
  ( main,
  )
where

import Control.Monad
import Mcmc
import Mcmc.Proposal
import Statistics.Distribution
import Statistics.Distribution.Normal
import System.Random.MWC hiding (uniform)

-- The state is composed of a tuple of two variables x and y.
type I = (Double, Double)

analysisName :: AnalysisName
analysisName = AnalysisName "Pair"

-- Improper prior for positive values.
pr :: I -> Log Double
pr (x, y) = largerThan 0.0001 x * largerThan 0.0001 y

-- Likelihood function.
lh :: I -> Log Double
lh (x, y) = exponential 1.0 (x + y)

-- Initial value.
start :: I
start = (1.1, 1.1)

-- Additive proposals with a custom Jacobian.
psX :: TuningParameter -> ProposalSimple I
psX t (x, y) g = do
  let d = normalDistr 0 t
  u <- genContinuous d g
  let x' = x + u
      qXY = Exp $ logDensity d u
      qYX = Exp $ logDensity d (negate u)
      r = abs $ (x + y) / (x + y + u)
      j = Exp $ log r
  return ((x', y), qYX / qXY, j)

psY :: TuningParameter -> ProposalSimple I
psY t (x, y) g = do
  let d = normalDistr 0 t
  u <- genContinuous d g
  let y' = y + u
      qXY = Exp $ logDensity d u
      qYX = Exp $ logDensity d (negate u)
      r = abs $ (x + y) / (x + y + u)
      j = Exp $ log r
  return ((x, y'), qYX / qXY, j)

prX :: Proposal I
prX =
  Proposal
    (PName "x")
    (PDescription "prX")
    (PDimension 1)
    (PWeight 1)
    (psX 1)
    (Just $ Tuner 1.0 psX)

prY :: Proposal I
prY =
  Proposal
    (PName "y")
    (PDescription "prY")
    (PDimension 1)
    (PWeight 1)
    (psY 1)
    (Just $ Tuner 1.0 psY)

-- -- The equivalent multiplicative proposals are commented out.

-- psX :: TuningParameter -> ProposalSimple I
-- psX t (x, y) g = do
--   let d = gammaDistr (recip t) t
--   u <- genContinuous d g
--   let x' = x * u
--       qXY = Exp $ logDensity d u
--       qYX = Exp $ logDensity d (recip u)
--       r = abs $ recip u * (x + y) / (x' + y)
--       j = Exp $ log r
--   return ((x', y), qYX / qXY, j)

-- psY :: TuningParameter -> ProposalSimple I
-- psY t (x, y) g = do
--   let d = gammaDistr (recip t) t
--   u <- genContinuous d g
--   let y' = y * u
--       qXY = Exp $ logDensity d u
--       qYX = Exp $ logDensity d (recip u)
--       r = abs $ recip u * (x + y) / (x + y')
--       j = Exp $ log r
--   return ((x, y'), qYX / qXY, j)

-- prX :: Proposal I
-- prX =
--   Proposal
--     (PName "x")
--     (PDescription "prX mult")
--     (PDimension 1)
--     (PWeight 1)
--     (psX 1)
--     (Just $ Tuner 1.0 psX)

-- prY :: Proposal I
-- prY =
--   Proposal
--     (PName "y")
--     (PDescription "prY mult")
--     (PDimension 1)
--     (PWeight 1)
--     (psY 1)
--     (Just $ Tuner 1.0 psY)

-- The proposal cycle consists of one proposal only. A uniform distribution is used to
-- slide the precision of the archer.
cc :: Cycle I
cc =
  cycleFromList
    [ prY
    , prX
    , slideContrarily 0 1.0 (PName "x y") (PWeight 1) Tune
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
          Debug
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg pr lh cc mon TraceAuto start g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings 3 1
  -- a <- mc3 mc3S pr lh cc mon (1, 1) g
  void $ mcmc mcmcS a
