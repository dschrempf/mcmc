{-# LANGUAGE OverloadedStrings #-}

{- |
Module      :  Archery
Description :  Modeling shots of an archer shots on a circular target
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May 27 17:11:28 2020.

Please refer to the corresponding [RevBayes
tutorial](https://revbayes.github.io/tutorials/mcmc/archery.html).

-}

module Main
  ( main
  ) where

import Control.Monad
import qualified Data.Text.Lazy         as T
import qualified Data.Text.Lazy.Builder as T
import qualified Data.Text.Lazy.Builder.RealFloat as T
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Exponential
import Statistics.Distribution.Gamma
import System.Random.MWC

import Statistics.Mcmc

type I = Double

-- Number of arrows.
n :: Int
n = 100

-- True precision of archer.
muTrue :: Double
muTrue = 1.0

-- Parameter of prior distribution.
alpha :: Double
alpha = 1.0

-- Likelihood function.
meanDistribution :: I -> GammaDistribution
-- meanDistribution x = gammaDistr (fromIntegral n) (fromIntegral n / x)
meanDistribution x = gammaDistr (fromIntegral n) (x / fromIntegral n)

-- Simulated mean distance from center.
arrowMean :: GenIO -> IO I
arrowMean = genContVar (meanDistribution muTrue)

priorDistribution :: ExponentialDistribution
priorDistribution = exponential alpha

negInf :: Fractional a => a
negInf = -(1/0)

prior :: I -> Log Double
prior x | x <= 0    = Exp negInf
        | otherwise = Exp $ logDensity priorDistribution x

likelihood :: I -> I -> Log Double
likelihood mu x | x <= 0    = Exp negInf
                | otherwise = Exp $ logDensity (meanDistribution mu) x

posterior :: I -> I -> Log Double
posterior mu x = prior x * likelihood mu x

moveCycle :: Cycle I
moveCycle = fromList [slideUniformDouble "mu; slide uniform" 1 1.0 True]

monRealFloat :: RealFloat a => MonitorParameter a
monRealFloat = MonitorParameter "Mu" (T.toStrict . T.toLazyText . T.formatRealFloat T.Fixed (Just 4))

monStd :: MonitorSimple I
monStd = MonitorSimple [monRealFloat] 1000

monFile :: MonitorFile I
monFile = monitorFile "Archery.log" [monRealFloat] 5

mons :: Monitor I
mons = Monitor monStd [monFile]

nBurn :: Maybe Int
nBurn = Just 2000

nAutoTune :: Maybe Int
nAutoTune = Just 400

nIter :: Int
nIter = 20000

main :: IO ()
main = do
  g <- create
  mu_observed <- arrowMean g
  putStrLn $ "True parameter: " <> show mu_observed
  let chain = mcmc 0 (posterior mu_observed) moveCycle mons g
  void $ mh nBurn nAutoTune nIter chain
