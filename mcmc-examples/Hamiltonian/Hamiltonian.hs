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

import Control.Lens
import Control.Monad
import qualified Data.Vector as VB
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Storable as VS
import Mcmc
import Numeric.AD (grad)
import qualified Numeric.LinearAlgebra as L
import System.Random.MWC hiding (uniform)

type IB = VB.Vector Double

type I = VS.Vector Double

-- instance ToJSON a => ToJSON (ZipList a)

-- instance FromJSON a => FromJSON (ZipList a)

-- The Rosenbrock function is really hard to estimate. It works somewhat, but
-- the example is a little bit frustrating.

-- logRosenbrock :: RealFloat a => a -> a -> [a] -> a
-- logRosenbrock a b [x, y] = log $ (a - x) ** 2 + b * ((y - x * x) ** 2)
-- logRosenbrock _ _ _ = error "lhf: Number of parameters has changed."

logGauss1 :: RealFloat a => a -> a -> a
logGauss1 s x = log (recip $ s * sqrt (2 * pi)) - 0.5 * x * x / s / s

logGaussN :: (VG.Vector v a, RealFloat a) => v a -> v a -> a
logGaussN ss xs = VG.foldl' (+) 0.0 (VG.zipWith logGauss1 ss xs)

standardDeviations :: (VG.Vector v a, Enum a, RealFloat a) => v a
-- Hard; dimension = 100.
standardDeviations = VG.fromList $ [0.02, 0.04 .. 1.0] ++ [2, 4 .. 100]
-- -- Easy; dimension = 10.
-- standardDeviations = VG.fromList [1, 2 .. 10]

dimension :: Int
dimension = VS.length (standardDeviations :: I)

llhf :: (VG.Vector v a, RealFloat a, Enum a) => v a -> a
llhf = logGaussN standardDeviations

gradientG :: IB -> IB
gradientG = grad llhf

gradient :: I -> I
gradient = VS.convert . gradientG . VS.convert

prf :: PriorFunction I
prf _ = 1.0

-- prf (ZipList [x, y]) = uniform (-2) 2 x * uniform (-1) 3 y
-- prf _ = error "prf: Number of parameters has changed."

-- The Likelihood function on I is less general, and does not allow for
-- automatic differentiation.
lhf :: LikelihoodFunction I
lhf = Exp . llhf

masses :: Masses
masses = L.trustSym $ L.diag $ L.fromList $ replicate dimension 1.0

hSettings :: HSettings I
hSettings =
  HSettings
    VB.convert
    (const VB.convert)
    gradient
    Nothing
    masses
    10
    0.05
    (HTune HTuneLeapfrog HTuneAllMasses)

initialState :: I
initialState = VS.fromList $ replicate dimension 1

hamiltonianProposal :: Proposal I
hamiltonianProposal = hamiltonian initialState hSettings n w
  where
    n = PName "Space"
    w = pWeight 1

cc :: Cycle I
cc =
  cycleFromList [hamiltonianProposal]

monPs :: [MonitorParameter I]
monPs = [view (singular (ix i)) >$< monitorDouble (n i) | i <- [0 .. (dimension -1)]]
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
          (BurnInWithCustomAutoTuning ([10, 20 .. 200] ++ [500, 500, 500]))
          (Iterations 10000)
          TraceAuto
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Info
  a <- mhg s prf lhf cc mon initialState g
  void $ mcmc s a
