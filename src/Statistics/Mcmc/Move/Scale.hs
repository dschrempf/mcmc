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

-- | Multiplicative move with Gamma distributed density.
scale
  :: Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> Double             -- ^ Shape.
  -> Double             -- ^ Scale.
  -> Bool               -- ^ Enable tuning.
  -> Move a
scale l n k th True  = moveGenericContinuous l n (gammaDistr k th) (*) (/)
                       (Just 1.0) $ Just (\t -> scale l n k (t*th) True)
scale l n k th False = moveGenericContinuous l n (gammaDistr k th) (*) (/)
                       Nothing Nothing

-- | Multiplicative move with Gamma distributed density; specialized to a one
-- dimensional state space of type 'Double'.
scaleDouble
  :: String             -- ^ Name.
  -> Double             -- ^ Shape.
  -> Double             -- ^ Scale.
  -> Bool               -- ^ Enable tuning.
  -> Move Double
scaleDouble = scale id

-- | Multiplicative move with Gamma distributed density. The scale of the Gamma
-- distributions is set to (shape)^{-1}, so that the mean of the Gamma
-- distribution is 1.0.
scaleUnbiased
  :: Lens' a Double     -- ^ Instruction about which parameter to change.
  -> String             -- ^ Name.
  -> Double             -- ^ Shape.
  -> Bool               -- ^ Enable tuning.
  -> Move a
scaleUnbiased l n k True  = moveGenericContinuous l n (gammaDistr k (1/k)) (*) (/)
                            (Just 1.0) $ Just (\t -> scaleUnbiased l n (k/t) True)
scaleUnbiased l n k False = moveGenericContinuous l n (gammaDistr k (1/k)) (*) (/)
                            Nothing Nothing
