{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

-- |
-- Module      :  Main
-- Description :  Simple tests
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- This example module is more involved and contains proposals with custom
-- Jacobians.
--
-- We sample the posterior distribution of a variable \(z \sim \mbox{Exp}(1.0)\)
-- separated into a sum \(z=x+y\), with \(x,y > 0\). The sampler propose new
-- values for \(x\) and \(y\).
module Main
  ( main,
  )
where

import Control.Lens
import Control.Monad
import Data.Aeson
import qualified Data.Vector.Fixed as VB
import qualified Data.Vector.Fixed.Boxed as VB
-- import qualified Data.Vector.Storable as VS
import Mcmc
-- import Numeric.AD.Double
-- import qualified Numeric.LinearAlgebra as L
import System.Random.MWC hiding (uniform)

-- The state is composed of two variables x and y.
type IG a = VB.Vec2 a

instance ToJSON a => ToJSON (VB.Vec2 a) where
  toJSON xs = toJSON [xs VB.! 0, xs VB.! 1]

type I = IG Double

analysisName :: AnalysisName
analysisName = AnalysisName "pair"

-- Improper prior for positive values of x and y. We set a bound strictly
-- greater than 0 because otherwise, numerical problems occur when the chain
-- traverses values very close to zero.
pr :: RealFloat a => PriorFunctionG (IG a) a
pr xs = greaterThan 0.00001 (xs VB.! 0) * greaterThan 0.00001 (xs VB.! 1)

f :: Num a => IG a -> a
-- f (x, y)= 3*x + 5*y
f xs = (xs VB.! 0) + (xs VB.! 1)

-- Likelihood function only acts on the sum of x and y, so we will need a custom
-- Jacobian in our proposals.
lh :: RealFloat a => LikelihoodFunctionG (IG a) a
lh = exponential 1.0 . f

-- Initial value.
start :: I
start = VB.mk2 1.1 1.1

-- Jacobian function. A lengthy derivation is required to find this function.
jacob :: RealFloat a => IG a -> Log a
jacob = Exp . log . recip . f

tupleL :: Lens' (VB.Vec2 a) (a, a)
tupleL = lens (\x -> (x VB.! 0, x VB.! 1)) (\_ (y1', y2') -> VB.mk2 y1' y2')

cc :: Cycle I
cc =
  cycleFromList
    [ -- The proposals require a custom Jacobian.
      liftProposalWith jacob (VB.element 0) (scaleUnbiased 1.0 (PName "x") (pWeight 1) Tune),
      liftProposalWith jacob (VB.element 1) (scaleUnbiased 1.0 (PName "y") (pWeight 1) Tune),
      -- Sliding the pair contrarily will not change the sum z, and so, no
      -- Jacobian is required (or the Jacobian will be 1.0). If
      -- 'scaleContrarily' is used or 'f' is not 1.0 for a contrary slide, the
      -- custom Jacobian is required.
      tupleL @~ slideContrarily 0 0.5 (PName "x y") (pWeight 3) Tune
    ]

-- TODO (high): Hamiltonian proposal with non-unit Jacobian. See note in the
-- Hamiltonian proposal module.

-- tspec :: HTuningSpec
-- tspec = either error id $ hTuningSpec masses 20 0.01 tconf
--   where
--     masses = L.trustSym $ L.ident $ L.size $ toVec start
--     tconf = HTuningConf HTuneLeapfrog HTuneAllMasses

-- toVec :: I -> VS.Vector Double
-- toVec xs = VS.generate 2 (\i -> xs VB.! i)

-- fromVec :: I -> VS.Vector Double -> I
-- fromVec _ xs = VB.mk2 (xs VS.! 0) (xs VS.! 1)

-- post :: RealFloat a => PosteriorFunctionG (IG a) a
-- post x = pr x * lh x

-- gradLogPosterior :: IG Double -> IG Double
-- gradLogPosterior = grad (ln . post)

-- hspec :: HSpec I
-- hspec = HSpec start toVec fromVec gradLogPosterior Nothing

-- hmc :: Proposal I
-- hmc = hamiltonian tspec hspec (PName "Hamiltonian") (pWeight 1)

-- cc :: Cycle I
-- cc = cycleFromList [hmc]

monPs :: [MonitorParameter I]
monPs =
  [ (VB.! 0) >$< monitorDouble "x",
    (VB.! 1) >$< monitorDouble "y",
    VB.sum >$< monitorDouble "x+y",
    VB.foldl1 (*) >$< monitorDouble "x*y",
    f >$< monitorDouble "f(x,y)"
  ]

monStd :: MonitorStdOut I
monStd = monitorStdOut (take 4 monPs) 1000

monFile :: MonitorFile I
monFile = monitorFile "" monPs 100

mon :: Monitor I
mon = Monitor monStd [monFile] []

main :: IO ()
main = do
  g <- create
  let mcmcS =
        Settings
          analysisName
          (BurnInWithAutoTuning 50000 500)
          (Iterations 200000)
          TraceAuto
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg mcmcS pr lh cc mon start g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings (NChains 4) (SwapPeriod 1) (NSwaps 1)
  -- a <- mc3 mc3S mcmcS pr lh cc mon TraceAuto start g
  void $ mcmc mcmcS a
