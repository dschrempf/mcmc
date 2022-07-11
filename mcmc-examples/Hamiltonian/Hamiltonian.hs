{-# OPTIONS_GHC -Wno-orphans #-}

-- |
-- Module      :  Main
-- Description :  Hamiltonian MCMC sampler
-- Copyright   :  2021 Dominik Schrempf
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
import qualified Data.Vector.Storable as VS
import Mcmc
import System.Random.Stateful

type IG = VB.Vector

type I = IG Double

-- -- The Rosenbrock function is really hard to estimate. It works somewhat, but
-- -- the example is a little bit frustrating.

-- logRosenbrock :: RealFloat a => a -> a -> IG a -> a
-- logRosenbrock a b xs
--   | n == 2 = log $ (a - x) ** 2 + b * ((y - x * x) ** 2)
--   | otherwise = error "lhf: Number of parameters has changed."
--   where
--     n = VB.length xs
--     x = xs VB.! 0
--     y = xs VB.! 1

-- lhf :: RealFloat a => LikelihoodFunctionG (IG a) a
-- lhf = Exp . negate . logRosenbrock 5 0.05

-- dimension :: Int
-- dimension = 2

logGauss1 :: RealFloat a => a -> a -> a
logGauss1 s x = log (recip $ s * sqrt (2 * pi)) - 0.5 * x * x / s / s

logGaussN :: RealFloat a => IG a -> IG a -> a
logGaussN ss xs = VB.foldl' (+) 0.0 (VB.zipWith logGauss1 ss xs)

standardDeviations :: RealFloat a => IG a
standardDeviations = VB.fromList $ map realToFrac xs
  where
    xs :: [Double]
    xs = [0.02, 0.04 .. 1.0] ++ [2, 4 .. 100]

lhf :: RealFloat a => LikelihoodFunctionG (IG a) a
lhf = Exp . logGaussN standardDeviations

dimension :: Int
dimension = VB.length (standardDeviations :: I)

initialState :: I
initialState = VB.fromList $ replicate dimension 1

htconf :: HTuningConf
htconf = HTuningConf HTuneLeapfrog HTuneAllMasses

hstruct :: HStructure IG
hstruct = HStructure initialState VS.convert (const VS.convert)

hTarget :: HTarget IG
hTarget = HTarget Nothing lhf Nothing

nutsProposal :: Proposal I
nutsProposal = nuts defaultNParams htconf hstruct hTarget n w
  where
    n = PName "Nuts"
    w = pWeight 1

cc :: Cycle I
cc =
  cycleFromList [nutsProposal]

monPs :: [MonitorParameter I]
monPs = [view (singular (ix i)) >$< monitorDouble (n i) | i <- [0 .. (dimension - 1)]]
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
  let g = mkStdGen 0
  let s =
        Settings
          (AnalysisName "hamiltonian")
          (BurnInWithCustomAutoTuning [] ([10, 20 .. 200] ++ replicate 20 200))
          (Iterations 8000)
          TraceAuto
          Overwrite
          Sequential
          Save
          LogStdOutAndFile
          Debug
  a <- mhg s noPrior lhf cc mon initialState g
  void $ mcmc s a
