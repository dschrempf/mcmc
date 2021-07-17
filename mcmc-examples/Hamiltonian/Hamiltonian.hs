{-# OPTIONS_GHC -Wno-orphans #-}

-- |
-- Module      :  Main
-- Description :  Hamiltonian MCMC sampler
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Mon Jul  5 17:11:32 2021.
module Main
  ( main,
  )
where

import Control.Applicative
import Control.Lens
import Control.Monad
import Data.Aeson
import Data.Foldable
import Mcmc
import Numeric.AD (grad)
import System.Random.MWC hiding (uniform)

type IG a = ZipList a

type I = IG Double

instance ToJSON a => ToJSON (ZipList a)

instance FromJSON a => FromJSON (ZipList a)

-- The Rosenbrock function is really hard to estimate. It works somewhat, but
-- the example is a little bit frustrating.

-- logRosenbrock :: RealFloat a => a -> a -> [a] -> a
-- logRosenbrock a b [x, y] = log $ (a - x) ** 2 + b * ((y - x * x) ** 2)
-- logRosenbrock _ _ _ = error "lhf: Number of parameters has changed."

logGauss1 :: RealFloat a => a -> a -> a
logGauss1 s x = log (recip $ s * sqrt (2 * pi)) - 0.5 * x * x / s / s

logGaussN :: RealFloat a => [a] -> [a] -> a
logGaussN ss xs = foldl' (+) 0.0 (zipWith logGauss1 ss xs)

standardDeviations :: (Enum a, RealFloat a) => [a]
standardDeviations = [0.02, 0.04 .. 1.0] ++ [2, 4 .. 100]

dimension :: Int
dimension = length (standardDeviations :: [Double])

llhf :: (RealFloat a, Enum a) => [a] -> a
llhf = logGaussN standardDeviations

gradient :: [Double] -> [Double]
gradient = grad llhf

prf :: PriorFunction I
prf _ = 1.0
-- prf (ZipList [x, y]) = uniform (-2) 2 x * uniform (-1) 3 y
-- prf _ = error "prf: Number of parameters has changed."

-- The Likelihood function on I is less general, and does not allow for
-- automatic differentiation.
lhfI :: LikelihoodFunction I
lhfI = Exp . llhf . getZipList

gradientI :: I -> I
gradientI = ZipList . gradient . getZipList

masses :: IG (Maybe Double)
masses = ZipList $ replicate dimension (Just 1)

hmcSettings :: HmcSettings ZipList
hmcSettings = HmcSettings gradientI masses 10 0.1 HmcTuneMassesAndLeapfrog

initialState :: I
initialState = ZipList $ replicate dimension 1

hmcProposal :: Proposal I
hmcProposal = hmc initialState hmcSettings n w
  where
    n = PName "Space"
    w = pWeight 1

cc :: Cycle I
cc =
  cycleFromList [hmcProposal]

monPs :: [MonitorParameter I]
monPs = [(view (singular (ix i)) . getZipList) >$< monitorDouble (n i) | i <- [0 .. (dimension -1)]]
  where
    n j = show j

mon :: Monitor I
mon =
  Monitor
    (monitorStdOut (take 4 monPs) 10)
    [monitorFile "params" monPs 2]
    []

main :: IO ()
main = do
  putStrLn $ "The dimension is: " <> show dimension <> "."
  g <- create
  let s =
        Settings
          (AnalysisName "hamiltonian")
          (BurnInWithCustomAutoTuning ([10, 20 .. 250] ++ [500, 500, 500]))
          (Iterations 10000)
          TraceAuto
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Debug
  a <- mhg s prf lhfI cc mon initialState g
  void $ mcmc s a
