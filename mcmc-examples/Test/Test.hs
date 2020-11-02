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
import Numeric.Log
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
proposals :: Cycle I
proposals =
  fromList
    [ _1 @~ slideSymmetric 1 "x" (Weight 5) Tune,
      _2 @~ slideSymmetric 1 "y" (Weight 5) Tune,
      scaleContrarily 1.0 1.0 "x y" (Weight 5) Tune
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

-- Number of burn in iterations.
nBurnIn :: Maybe Int
nBurnIn = Just 20000

-- Auto tuning period.
nAutoTune :: Maybe Int
nAutoTune = Just 1000

-- Number of Metropolis-Hastings iterations after burn in.
nIter :: Int
nIter = 100000

main :: IO ()
main = do
  g <- create
  let s = force $ status "test" pr lh proposals mon (1.0, 1.0) nBurnIn nAutoTune nIter g
  void $ mh s
