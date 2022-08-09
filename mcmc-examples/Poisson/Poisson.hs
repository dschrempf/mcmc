-- |
-- Module      :  Main
-- Description :  Estimate airline fatal accidents
-- Copyright   :  2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Mon Aug  8 15:19:33 2022.
--
-- See my [blog post about estimating airline fatal
-- accidents](http://dschrempf.github.io/coding/2022-06-28-sample-from-a-posterior-using-markov-chain-monte-carlo-algorithms-and-haskell/).
module Main
  ( main,
  )
where

-- We need 'void'.
import Control.Monad
-- I am developing the 'mcmc' library.
import Mcmc
-- We need to sample random numbers; requires the 'random' library.
import System.Random

-- The state space 'I' is the fatal accident rate 'lambda', a floating point
-- number.
type I = Double

lh :: LikelihoodFunction I
lh x = product [poisson x y | y <- fatalities]
  where
    fatalities = [24, 25, 31, 31, 22, 21, 26, 20, 16, 22]

ppSl :: Proposal I
ppSl = slideSymmetric 0.1 (PName "lambda") (pWeight 1) Tune

ppSc :: Proposal I
ppSc = scaleUnbiased 0.1 (PName "lambda") (pWeight 1) Tune

cc :: Cycle I
cc = cycleFromList [ppSl, ppSc]

-- 'monitorDouble' is a simple monitor printing the value of a 'Double'.
monLambda :: MonitorParameter I
monLambda = monitorDouble "lambda"

-- We print the value of lambda to the standard output every 100 iterations.
monStdOut :: MonitorStdOut I
monStdOut = monitorStdOut [monLambda] 100

-- We log the value of lambda to a file more often.
monFile :: MonitorFile I
monFile = monitorFile "lambda" [monLambda] 3

mon :: Monitor I
mon = Monitor monStdOut [monFile] []

ss :: Settings
ss =
  Settings
    -- Provide an analysis name.
    (AnalysisName "poisson")
    -- Burn in for 2000 generations. During burn in, the proposals are tuned
    -- automatically. This is called "auto tuning". Here, auto tuning is
    -- performed every 200 iterations.
    (BurnInWithAutoTuning 2000 200)
    -- Number of actual iterations after burn in.
    (Iterations 30000)
    -- The trace of the Markov chain contains the attained values. In our case,
    -- it is a vector of fatal accident rates. Here, we tell the sampler to use
    -- the shortest trace possible. In our case, this will be a single value.
    -- However, when using batch monitors, or when auto tuning the masses of
    -- proposals based on Hamiltonian dynamics, the required length of the trace
    -- is larger than 1. masses. The trace length can also be set manually.
    TraceAuto
    -- Overwrite files created by a possible previous analysis.
    Overwrite
    -- Do not run chains in parallel. For the standard Metropolis-Hastings-Green
    -- algorithm, this has no effect. However, there are algorithms such as the
    -- MC3 algorithm with multiple chains that can run in parallel.
    Sequential
    -- Save the chain so that it can be continued (see 'mcmcContinue').
    Save
    -- Log to standard output and save the log to a file.
    LogStdOutAndFile
    -- Verbosity.
    Info

main :: IO ()
main = do
  let g = mkStdGen 0
  -- Set up the Markov chain. For computational efficiency (mutable vectors),
  -- this requires IO.
  al <- mhg ss noPrior lh cc mon 1.0 g
  -- We ignore the actual return value which is the complete Markov chain object
  -- using 'void'.
  void $ mcmc ss al
