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
import Control.Monad
import Data.Aeson
import Mcmc
import Numeric.AD (grad)
import System.Random.MWC hiding (uniform)

type I = ZipList Double

type IL = ZipList (Log Double)

instance ToJSON a => ToJSON (ZipList a)

instance FromJSON a => FromJSON (ZipList a)

rosenbrock :: RealFloat a => a -> a -> [a] -> a
rosenbrock a b [x, y] = log $ (a - x) ** 2 + b * ((y - x * x) ** 2)
rosenbrock _ _ _ = error "lhf: Number of parameters has changed."

rb :: RealFloat a => [a] -> a
rb = rosenbrock 1 100

gradient :: [Double] -> [Log Double]
gradient = map (Exp . log) . grad rb

prf :: PriorFunction I
prf (ZipList [x, y]) = uniform (-2) 2 x * uniform (-1) 3 y
prf _ = error "prf: Number of parameters has changed."

-- The Llkelihood function on I is less general, and does not allow for
-- automatic differentiation.
lhfI :: LikelihoodFunction I
lhfI = Exp . rb . getZipList

gradientI :: I -> IL
gradientI = ZipList . gradient . getZipList

masses :: I
masses = ZipList [1, 1]

hmcSettings :: HmcSettings ZipList
hmcSettings = HmcSettings gradientI masses 10 0.1

initialState :: I
initialState = ZipList [1.1, 0.9]

hmcProposal :: Proposal I
hmcProposal = hmc initialState hmcSettings n w Tune
  where
    n = PName "Space"
    w = pWeight 1

cc :: Cycle I
cc =
  cycleFromList [hmcProposal]

monPs :: [MonitorParameter I]
monPs = [head . getZipList >$< monitorDouble "x", head . tail . getZipList >$< monitorDouble "y"]

mon :: Monitor I
mon =
  Monitor
    (monitorStdOut monPs 10)
    [monitorFile "params" monPs 2]
    []

main :: IO ()
main = do
  g <- create
  let s =
        Settings
          (AnalysisName "hamiltonian")
          (BurnInWithCustomAutoTuning [10,20..250])
          (Iterations 10000)
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Debug
  a <- mhg prf lhfI cc mon TraceAuto initialState g
  void $ mcmc s a
