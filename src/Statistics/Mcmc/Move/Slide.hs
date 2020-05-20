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
import Statistics.Mcmc.Move

-- | Additive move with normally distributed density.
slide
  :: Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> Double             -- ^ Mean.
  -> Double             -- ^ Standard deviation.
  -> Bool               -- ^ Enable tuning.
  -> Move a
slide l n m s True  = moveGenericContinuous l n (normalDistr m s) (+) (-)
                      (Just 1.0) $ Just (\t -> slide l n m (t*s) True)
slide l n m s False = moveGenericContinuous l n (normalDistr m s) (+) (-)
                      Nothing Nothing

-- | Additive move with normally distributed density; specialized to a one
-- dimensional state space of type 'Double'.
slideDouble
  :: String -- ^ Name.
  -> Double -- ^ Mean.
  -> Double -- ^ Standard deviation.
  -> Bool   -- ^ Enable tuning.
  -> Move Double
slideDouble = slide id
