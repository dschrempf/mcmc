{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Statistics.Mcmc.Move.Slide
Description :  Normally distributed move
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 10:59:13 2020.

-}

module Statistics.Mcmc.Move.Slide
  ( slide
  , slideDouble
  ) where

import Lens.Micro
import Statistics.Distribution.Normal

import Statistics.Mcmc.Move.Generic
import Statistics.Mcmc.Move.Types

-- | Additive move with normally distributed density.
slide
  :: Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> Double             -- ^ Mean.
  -> Double             -- ^ Standard deviation.
  -> Move a
slide l n m s = moveGenericContinuous l n (normalDistr m s) (+) (-)

-- | Additive move with normally distributed density; specialized to a one
-- dimensional state space of type 'Double'.
slideDouble :: String -> Double -> Double -> Move Double
slideDouble = slide id
