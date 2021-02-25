-- |
-- Module      :  Poisson
-- Description :  Poisson regression model for airline fatalities
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May  8 13:30:51 2020.
--
-- See https://revbayes.github.io/tutorials/mcmc/poisson.html.
module Poisson
  ( poissonBench,
  )
where

import Control.Monad
import Lens.Micro
import Mcmc
import System.Random.MWC

type I = (Double, Double)

fatalities :: [Int]
fatalities = [24, 25, 31, 31, 22, 21, 26, 20, 16, 22]

normalizedYears :: [Double]
normalizedYears = map (subtract m) ys
  where
    ys = [1976.0 .. 1985.0]
    m = sum ys / fromIntegral (length ys)

f :: Int -> Double -> I -> Log Double
f ft yr (a, b) = poisson l (fromIntegral ft)
  where
    l = exp $ a + b * yr

lh :: LikelihoodFunction I
lh x = product [f ft yr x | (ft, yr) <- zip fatalities normalizedYears]

proposalAlpha :: Proposal I
proposalAlpha = _1 @~ slideSymmetric 0.2 (PName "Alpha") (PWeight 1) NoTune

proposalBeta :: Proposal I
proposalBeta = _2 @~ slideSymmetric 0.2 (PName "Beta") (PWeight 1) NoTune

proposals :: Cycle I
proposals = cycleFromList [proposalAlpha, proposalBeta]

initial :: I
initial = (0, 0)

monAlpha :: MonitorParameter I
monAlpha = fst >$< monitorDouble "alpha"

monBeta :: MonitorParameter I
monBeta = snd >$< monitorDouble "beta"

monStd :: MonitorStdOut I
monStd = monitorStdOut [monAlpha, monBeta] 150

mon :: Monitor I
mon = Monitor monStd [] []

poissonBench :: GenIO -> IO ()
poissonBench g = do
  let s =
        Settings
          (AnalysisName "Poisson")
          (BurnInWithAutoTuning 2000 200)
          (Iterations 10000)
          Overwrite
          Sequential
          NoSave
          LogStdOutOnly
          Quiet
  a <- mhg noPrior lh proposals mon TraceAuto initial g
  void $ mcmc s a
