{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Main
-- Description :  Estimate the parameter vector of a Dirichlet distribution
-- Copyright   :  (c) Dominik Schrempf, 2020
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
import Statistics.Distribution.Dirichlet
import System.Random.MWC

data I = I
  { _alphas :: Simplex,
    _norm :: Double
  }
  deriving (Eq, Show, Generic)

-- Create accessors (lenses) to the parameters in the state space.Nothing
makeLenses ''I

instance ToJSON I

instance FromJSON I

-- The true parameter values.
alphasTrue :: V.Vector Double
alphasTrue = V.fromList [0.1, 10, 8.0, 4.2, 4.5, 0.3]

nObservations :: Int
nObservations = 100

-- Simulate data using a Dirichlet distribution with the true parameter values.
simulateData :: GenIO -> IO [V.Vector Double]
simulateData g = replicateM nObservations (dirichletSample dd g)
  where
    dd = either error id $ dirichletDistribution alphasTrue

-- Prior function.
pr :: PriorFunction I
pr (I as n)
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

alphaProposals :: [Proposal I]
alphaProposals =
  [ alphas @~ beta i (PName "Alpha") (PWeight 1) Tune
    | i <- [0 .. (V.length alphasTrue - 1)]
  ]

-- -- Cycle with Dirichlet proposal.
-- proposals :: Cycle I
-- proposals =
--   cycleFromList
--     [ alphas @~ dirichlet "dirichlet" 1 True,
--       norm @~ scaleUnbiased 8.0 "scale norm" 1 True
--     ]

-- Cycle with beta proposals.
cc :: Cycle I
cc =
  cycleFromList $
    norm @~ scaleUnbiased 8.0 (PName "Norm") (PWeight 1) Tune :
    alphaProposals

monNorm :: MonitorParameter I
monNorm = _norm >$< monitorDouble "Norm"

monN :: Int -> MonitorParameter I
monN n = (\x -> toVector (_alphas x) V.! n * _norm x) >$< monitorDouble name
  where
    name = "Alpha " <> show n

monAlphas :: [MonitorParameter I]
monAlphas = [monN i | i <- [0 .. (V.length alphasTrue - 1)]]

monStdOut :: MonitorStdOut I
monStdOut = monitorStdOut [monNorm, monN 0, monN 1, monN 2] 100

monFile :: MonitorFile I
monFile = monitorFile "" (monNorm : monAlphas) 50

mon :: Monitor I
mon = Monitor monStdOut [monFile] []

-- Starting state.
start :: I
start = I as 1.0
  where
    as = simplexUniform (V.length alphasTrue)

burnIn :: BurnIn
burnIn = BurnInWithAutoTuning 3000 100

iterations :: Int
iterations = 30000

main :: IO ()
main = do
  g <- create
  putStrLn "-- Simulate Data."
  xs <- simulateData g
  print xs
  print start
  let
    lh = lhf xs
    s = Settings "dirichlet" burnIn iterations Overwrite NoSave Info
    a = mhg pr lh cc mon start g
  _ <- mcmc s a
  putStrLn "Done."
