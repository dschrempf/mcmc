-- |
-- Module      :  Main
-- Description :  Gaussian model; test marginal likelihood
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Jan 12 13:00:40 2021.
--
-- See "A Gaussian Model" in Lartillot, N., & Philippe, H., Computing Bayes
-- Factors Using Thermodynamic Integration, Systematic Biology, 55(2), 195–207
-- (2006). http://dx.doi.org/10.1080/10635150500433722.
module Main
  ( main,
  )
where

import Control.Lens
import Mcmc
import System.Random.MWC

-- Variance.
v :: Double
v = 0.01

-- Dimension.
d :: Int
d = 1

-- Sample size per point.
n :: Int
n = 2000

-- State space.
type I = [Double]

-- Prior function.
prf :: PriorFunction I
prf = product . map (normal 0 1)

-- Gaussian exponential function.
gef :: Double -> Log Double
gef x = Exp $ negate $ x * x / (2 * v)

-- Likelihood function.
lhf :: LikelihoodFunction I
lhf = product . map gef

-- Cycle.
cc :: Cycle I
cc =
  cycleFromList
    [ singular (ix i)
        @~ slideUniformSymmetric 1.0 (PName $ "x" ++ show i) (pWeight 1) Tune
      | i <- [0 .. (d -1)]
    ]

main :: IO ()
main = do
  g <- create
  let ss =
        MLSettings
          (AnalysisName "gauss")
          SteppingStoneSampling
          (NPoints 512)
          (BurnInWithAutoTuning 5000 200)
          (BurnInWithAutoTuning 1000 100)
          (Iterations n)
          Overwrite
          LogStdOutAndFile
          Info
      i0 = replicate d 0
  _ <- marginalLikelihood ss prf lhf cc (simpleMonitor 1000) i0 g
  putStrLn $ "The correct value is: "  ++ show corVal
  where corVal = 0.5 * fromIntegral d * (log v - log (1+v))
