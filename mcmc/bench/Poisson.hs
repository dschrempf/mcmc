{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

-- |
-- Module      :  Poisson
-- Description :  Poisson regression model for airline fatalities
-- Copyright   :  2021 Dominik Schrempf
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
    poissonHamiltonianBench,
  )
where

import Control.Monad
import Data.Aeson
import Data.Typeable
import qualified Data.Vector.Fixed as VB
import qualified Data.Vector.Fixed.Boxed as VB
import qualified Data.Vector.Storable as VS
import Mcmc
import System.Random.Stateful

type IG = VB.Vec2

type I = IG Double

instance (VB.Arity n, ToJSON a) => ToJSON (VB.Vec n a) where
  toJSON = toJSON . VB.toList

fatalities :: [Int]
fatalities = [24, 25, 31, 31, 22, 21, 26, 20, 16, 22]

normalizedYears :: [Double]
normalizedYears = map (subtract m) ys
  where
    ys = [1976.0 .. 1985.0]
    m = sum ys / fromIntegral (length ys)

f :: (RealFloat a, Typeable a) => Int -> Double -> IG a -> LikelihoodG a
f ft yr xs = poisson l ft
  where
    a = xs VB.! 0
    b = xs VB.! 1
    l = exp $ a + b * realToFrac yr

lh :: (RealFloat a, Typeable a) => LikelihoodFunctionG (IG a) a
lh x = product [f ft yr x | (ft, yr) <- zip fatalities normalizedYears]

proposalAlpha :: Proposal I
proposalAlpha = (VB.element 0) @~ slideSymmetric 0.2 (PName "Alpha") (pWeight 1) NoTune

proposalBeta :: Proposal I
proposalBeta = (VB.element 1) @~ slideSymmetric 0.2 (PName "Beta") (pWeight 1) NoTune

proposals :: Cycle I
proposals = cycleFromList [proposalAlpha, proposalBeta]

initial :: I
initial = VB.mk2 0 0

monAlpha :: MonitorParameter I
monAlpha = (VB.! 0) >$< monitorDouble "alpha"

monBeta :: MonitorParameter I
monBeta = (VB.! 1) >$< monitorDouble "beta"

monStd :: MonitorStdOut I
monStd = monitorStdOut [monAlpha, monBeta] 150

-- monF :: MonitorFile I
-- monF = monitorFile "params" [monAlpha, monBeta] 10

mon :: Monitor I
mon = Monitor monStd [] []

poissonBench :: IOGenM StdGen -> IO ()
poissonBench g = do
  let s =
        Settings
          (AnalysisName "PoissonR")
          (BurnInWithAutoTuning 2000 200)
          (Iterations 10000)
          TraceAuto
          Overwrite
          Sequential
          NoSave
          LogStdOutOnly
          Quiet
  a <- mhg s noPrior lh proposals mon initial g
  void $ mcmc s a

toVec :: I -> VS.Vector Double
toVec xs = VS.generate 2 (\i -> xs VB.! i)

fromVec :: I -> VS.Vector Double -> I
fromVec _ xs = VB.mk2 (xs VS.! 0) (xs VS.! 1)

hparams :: HParams
hparams = HParams Nothing Nothing Nothing

htconf :: HTuningConf
htconf = HTuningConf HTuneLeapfrog HTuneAllMasses

hstruct :: HStructure IG
hstruct = HStructure initial toVec fromVec

target :: HTarget IG
target = HTarget Nothing lh Nothing

hmcProposal :: Cycle I
hmcProposal = cycleFromList [hamiltonian hparams htconf hstruct target (PName "Hamiltonian") (pWeight 1)]

poissonHamiltonianBench :: IOGenM StdGen -> IO ()
poissonHamiltonianBench g = do
  let s =
        Settings
          (AnalysisName "PoissonH")
          (BurnInWithAutoTuning 2000 200)
          (Iterations 10000)
          TraceAuto
          Overwrite
          Sequential
          NoSave
          LogStdOutOnly
          Quiet
  a <- mhg s noPrior lh hmcProposal mon initial g
  void $ mcmc s a
