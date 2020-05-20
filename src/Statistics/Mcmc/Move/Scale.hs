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
  , scaleNeutral
  ) where

import Lens.Micro
import Statistics.Distribution.Gamma

import Statistics.Mcmc.Move.Generic
import Statistics.Mcmc.Move.Types

-- | Multiplicative move with Gamma distributed density.
scale
  :: Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> Double             -- ^ Shape.
  -> Double             -- ^ Scale.
  -> Move a
scale l n k t = moveGenericContinuous l n (gammaDistr k t) (*) (/)

-- | Multiplicative move with Gamma distributed density; specialized to a one
-- dimensional state space of type 'Double'.
scaleDouble :: String -> Double -> Double -> Move Double
scaleDouble = scale id

-- | Multiplicative move with Gamma distributed density. The scale of the Gamma
-- distributions is set to (shape)^{-1}, so that the mean of the Gamma
-- distribution is 1.0.
--
-- XXX: "Neutral" is probably not the best name, can we think of a better
-- alternative?.
scaleNeutral
  :: Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> Double             -- ^ Shape.
  -> Move a
scaleNeutral l n k = moveGenericContinuous l n (gammaDistr k (1/k)) (*) (/)
