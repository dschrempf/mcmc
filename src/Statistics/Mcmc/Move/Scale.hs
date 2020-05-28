{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Statistics.Mcmc.Move.Scale
Description :  Scaling move with Gamma distribution
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 14 21:49:23 2020.

-}

module Statistics.Mcmc.Move.Scale
  ( scale
  , scaleDouble
  , scaleUnbiased
  ) where

import Lens.Micro
import Statistics.Distribution.Gamma

import Statistics.Mcmc.Move.Generic
import Statistics.Mcmc.Move

-- The actual move with tuning parameter.
scaleSimple
  :: Lens' a Double
  -> Double
  -> Double
  -> Double
  -> MoveSimple a
scaleSimple l k th t = moveGenericContinuous l (gammaDistr k (t*th)) (*) (/)

-- | Multiplicative move with Gamma distributed density.
scale
  :: String             -- ^ Name.
  -> Int                -- ^ Weight.
  -> Lens' a Double     -- ^ Instruction about which parameter to change.
  -> Double             -- ^ Shape.
  -> Double             -- ^ Scale.
  -> Bool               -- ^ Enable tuning.
  -> Move a
scale n w l k th True  = Move n w (scaleSimple l k th 1.0) (Just $ tuner $ scaleSimple l k th)
scale n w l k th False = Move n w (scaleSimple l k th 1.0)  Nothing

-- | Multiplicative move with Gamma distributed density; specialized to a one
-- dimensional state space of type 'Double'.
scaleDouble
  :: String             -- ^ Name.
  -> Int                -- ^ Weight.
  -> Double             -- ^ Shape.
  -> Double             -- ^ Scale.
  -> Bool               -- ^ Enable tuning.
  -> Move Double
scaleDouble n w = scale n w id

-- | Multiplicative move with Gamma distributed density. The scale of the Gamma
-- distributions is set to (shape)^{-1}, so that the mean of the Gamma
-- distribution is 1.0.
scaleUnbiased
  :: String             -- ^ Name.
  -> Int                -- ^ Weight.
  -> Lens' a Double     -- ^ Instruction about which parameter to change.
  -> Double             -- ^ Shape.
  -> Bool               -- ^ Enable tuning.
  -> Move a
scaleUnbiased n w l k True  = Move n w (scaleSimple l k (1/k) 1.0) (Just $ tuner $ scaleSimple l k (1/k))
scaleUnbiased n w l k False = Move n w (scaleSimple l k (1/k) 1.0)  Nothing
