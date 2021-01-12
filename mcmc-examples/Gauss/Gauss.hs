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

import Mcmc

-- Variance.
v :: Double
v = 1

-- Dimension.
d :: Int
d = 1

-- Sample size.
n :: Int
n = 1000

type I = [Double]

ef :: Double -> Double
ef x = exp (negate $ x * x / (2 * v))

lhf :: LikelihoodFunction I
lhf = product . map (Exp . log . ef)

main :: IO ()
main = undefined
