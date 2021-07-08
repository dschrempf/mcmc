{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Main
-- Description :  Estimate the parameter vector of a Dirichlet distribution
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Oct 20 14:04:04 2020.
module Main
  ( main,
  )
where

import Control.Lens
import Control.Monad
import Data.Aeson
import qualified Data.Vector.Unboxed as V
import GHC.Generics
import Mcmc
-- import Statistics.Distribution
import Statistics.Distribution.Dirichlet
-- import Statistics.Distribution.Normal
import System.Random.MWC

-- The Dirichlet distribution is parametrized by a set of alpha values.
alphasTrue :: V.Vector Double
alphasTrue = V.fromList $ map (/ 10) [1 .. 25]

-- The state space includes the scaled alpha parameters with sum 1.0, and their
-- normalization constant.
data I = I
  { _alphas :: Simplex,
    _norm :: Double
  }
  deriving (Eq, Show, Generic)

-- Create accessors (lenses) to the parameters in the state space.
makeLenses ''I

-- Make sure that we can store and restore the Markov chain.
instance ToJSON I

instance FromJSON I

nObservations :: Int
nObservations = 800

-- Simulate data using a Dirichlet distribution with the true parameter values.
simulateData :: GenIO -> IO [V.Vector Double]
simulateData g = replicateM nObservations (dirichletSample dd g)
  where
    dd = either error id $ dirichletDistribution alphasTrue

-- Improper, uninformative prior function. We don't allow negative values.
prf :: PriorFunction I
prf (I as n)
  | V.any (< 0) (toVector as) = 0.0
  | n < 0 = 0.0
  | otherwise = 1.0

-- The likelihood function is just the product of the Dirichlet densities of all
-- observations.
lhf :: [V.Vector Double] -> LikelihoodFunction I
lhf xs (I as n) = case eitherDds of
  Left _ -> 0
  Right dd -> product $ map (dirichletDensity dd) xs
  where
    eitherDds = dirichletDistribution $ V.map (* n) $ toVector as

-- Beta proposals on the simplex storing the alpha parameters.
alphaProposals :: [Proposal I]
alphaProposals =
  [ alphas @~ beta i (PName "Alpha") (pWeight 1) Tune
    | i <- [0 .. (V.length alphasTrue - 1)]
  ]

-- Cycle with the beta proposals and a proposal changing the normalization
-- constant.
cc :: Cycle I
cc =
  cycleFromList $
    norm @~ scaleUnbiased 8.0 (PName "Norm") (pWeight 1) Tune :
    alphaProposals

-- -- Cycle with Dirichlet proposal.
-- proposals :: Cycle I
-- proposals =
--   cycleFromList
--     [ alphas @~ dirichlet "dirichlet" 1 True,
--       norm @~ scaleUnbiased 8.0 "scale norm" 1 True
--     ]

monNorm :: MonitorParameter I
monNorm = _norm >$< monitorDouble "Norm"

monN :: Int -> MonitorParameter I
monN n = (\x -> toVector (_alphas x) V.! n * _norm x) >$< monitorDouble name
  where
    name = "Alpha " <> show n

monAlphas :: [MonitorParameter I]
monAlphas = [monN i | i <- [0 .. (V.length alphasTrue - 1)]]

-- Careful, this only works when the dimension is three or larger.
monStdOut :: MonitorStdOut I
monStdOut = monitorStdOut [monNorm, monN 0, monN 1, monN 2] 100

monFile :: MonitorFile I
monFile = monitorFile "" (monNorm : monAlphas) 2

mon :: Monitor I
mon = Monitor monStdOut [monFile] []

-- Starting state.
start :: I
start = I as s
  where
    -- Bad starting value.
    s = 1.0
    as = simplexUniform (V.length alphasTrue)
    -- -- Good starting value.
    -- s = V.sum alphasTrue
    -- normalize = V.map (/ s)
    -- as = either error id $ simplexFromVector $ normalize alphasTrue

main :: IO ()
main = do
  g <- create
  putStrLn "-- Simulate Data."
  xs <- simulateData g
  -- print xs
  -- print start
  let -- Settings.
      s =
        Settings
          (AnalysisName "dirichlet")
          (BurnInWithCustomAutoTuning [100, 100, 200, 300, 400, 500, 500, 500, 500])
          (Iterations 10000)
          TraceAuto
          Overwrite
          Parallel
          NoSave
          LogStdOutAndFile
          Info
  -- Initialize the Metropolis-Hastings-Green algorithm.
  a <- mhg s prf (lhf xs) cc mon start g
  -- -- Metropolis-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings (NChains 3) (SwapPeriod 2) (NSwaps 1)
  -- a <- mc3 mc3S prf (lhf xs) cc mon TraceAuto start g
  -- Run the MCMC sampler.
  _ <- mcmc s a
  putStrLn "Done."
