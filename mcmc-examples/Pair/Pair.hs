{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE RankNTypes #-}

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

-- import Control.Lens
import Control.Monad
import qualified Data.Vector as VB
import qualified Data.Vector.Storable as VS
import Mcmc
import qualified Numeric.LinearAlgebra as L
import System.Random.MWC hiding (uniform)

-- The state is composed of a vector containing two variables x and y.
type IG = VB.Vector

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
start = VB.fromList [1.1, 1.1]

-- Jacobian function. A lengthy derivation is required to find this function.
jacob :: RealFloat a => IG a -> Log a
jacob = Exp . log . recip . f

-- -- The lens functions could be improved. How?
-- tupleL :: Lens' (IG a) (a, a)
-- tupleL = lens (\x -> (x VB.! 0, x VB.! 1)) (\_ (y1', y2') -> VB.fromList [y1', y2'])

-- -- The lens functions could be improved. How?
-- indexL :: Int -> Lens' (IG a) a
-- indexL i = singular (ix i)

-- cc :: Cycle I
-- cc =
--   cycleFromList
--     [ -- The proposals require a custom Jacobian.
--       liftProposalWith jacob (indexL 0) (scaleUnbiased 1.0 (PName "x") (pWeight 1) Tune),
--       liftProposalWith jacob (indexL 1) (scaleUnbiased 1.0 (PName "y") (pWeight 1) Tune),
--       -- Sliding the pair contrarily will not change the sum z, and so, no
--       -- Jacobian is required (or the Jacobian will be 1.0). If
--       -- 'scaleContrarily' is used or 'f' is not 1.0 for a contrary slide, the
--       -- custom Jacobian is required.
--       tupleL @~ slideContrarily 0 0.5 (PName "x y") (pWeight 3) Tune
--     ]

tspec :: HTuningSpec
tspec = either error id $ hTuningSpec masses 10 0.1 tconf
  where
    masses = L.trustSym $ L.scale 5 $ L.matrix 2 [1.0, 0.0, 0.0, 1.0]
    tconf = HTuningConf HNoTuneLeapfrog HTuneDiagonalMassesOnly

toVec :: I -> VS.Vector Double
toVec = VS.convert

fromVec :: I -> VS.Vector Double -> I
fromVec = const VS.convert

hspec :: HSpec IG
hspec = HSpec start toVec fromVec

htarget :: HTarget IG
htarget = HTarget (Just pr) lh (Just jacob)

hmc :: Proposal I
hmc = hamiltonian tspec hspec htarget (PName "Hamiltonian") (pWeight 1)

cc :: Cycle I
cc = cycleFromList [hmc]

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
monFile = monitorFile "" monPs 1

mon :: Monitor I
mon = Monitor monStd [monFile] []

main :: IO ()
main = do
  g <- create
  let mcmcS =
        Settings
          analysisName
          (BurnInWithCustomAutoTuning [] $ [10, 20 .. 200] ++ replicate 10 500)
          (Iterations 20000)
          TraceAuto
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Debug
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg mcmcS pr lh cc mon start g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings (NChains 4) (SwapPeriod 1) (NSwaps 1)
  -- a <- mc3 mc3S mcmcS pr lh cc mon TraceAuto start g
  void $ mcmc mcmcS a
