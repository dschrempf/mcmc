-- |
-- Module      :  Main
-- Description :  Simple tests
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- This example module is more involved and contains proposals with custom
-- Jacobians.
--
-- We sample the posterior distribution of a variable \(z \sim \mbox{Exp}(1.0)\)
-- separated into a sum \(z=x+y\), with \(x,y > 0\). The sampler propose new
-- values for \(x\) and \(y\).
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
analysisName = AnalysisName "pair"

-- Improper prior for positive values of x and y. We set a bound strictly
-- greater than 0 because otherwise, numerical problems occur when the chain
-- traverses values very close to zero.
pr :: PriorFunction I
pr (x, y) = greaterThan 0.00001 x * greaterThan 0.00001 y

f :: I -> Double
-- f (x, y)= 3*x + 5*y
f (x, y)= x + y

-- Likelihood function only acts on the sum of x and y, so we will need a custom
-- Jacobian in our proposals.
lh :: LikelihoodFunction I
lh = exponential 1.0 . f

-- Initial value.
start :: I
start = (1.1, 1.1)

-- Jacobian function. A lengthy derivation is required to find this function.
jacobian :: JacobianFunction I
jacobian = Exp . log . recip . f

cc :: Cycle I
cc =
  cycleFromList
    [
      -- The proposals require a custom Jacobian.
      liftProposalWith jacobian _1 (scaleUnbiased 1.0 (PName "x") (pWeight 1) Tune),
      liftProposalWith jacobian _2 (scaleUnbiased 1.0 (PName "y") (pWeight 1) Tune),
      -- Sliding the pair contrarily will not change the sum z, and so, no
      -- Jacobian is required (or the Jacobian will be 1.0). If
      -- 'scaleContrarily' is used or 'f' is not 1.0 for a contrary slide, the
      -- custom Jacobian is required.
      slideContrarily 0 0.5 (PName "x y") (pWeight 3) Tune
    ]

monPs :: [MonitorParameter I]
monPs =
  [ fst >$< monitorDouble "x",
    snd >$< monitorDouble "y",
    uncurry (+) >$< monitorDouble "x+y",
    uncurry (*) >$< monitorDouble "x*y",
    f >$< monitorDouble "f(x,y)"
  ]

monStd :: MonitorStdOut I
monStd = monitorStdOut (take 4 monPs) 1000

monFile :: MonitorFile I
monFile = monitorFile "" monPs 100

mon :: Monitor I
mon = Monitor monStd [monFile] []

main :: IO ()
main = do
  g <- create
  let mcmcS =
        Settings
          analysisName
          (BurnInWithAutoTuning 50000 500)
          (Iterations 200000)
          TraceAuto
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg mcmcS pr lh cc mon start g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings (NChains 4) (SwapPeriod 1) (NSwaps 1)
  -- a <- mc3 mc3S mcmcS pr lh cc mon TraceAuto start g
  void $ mcmc mcmcS a
