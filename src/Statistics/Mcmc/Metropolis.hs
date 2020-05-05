{- |
Module      :  Statistics.Mcmc.Metropolis
Description :  Metropolis-Hastings at its best
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 20:11:30 2020.

-}

module Statistics.Mcmc.Metropolis
  (
  ) where

import Control.Monad
import Control.Monad.Primitive
import Numeric.Log
import System.Random.MWC

import Statistics.Mcmc.Types

{-# INLINE mhRatio #-}
mhRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
mhRatio lX lY qXY qYX = lY + qYX - lX - qXY

mhMove :: PrimMonad m
  -- ^ Current state and log-likelihood as well as the log-likelihood function
  -- and a generator.
  => (Link a, LogLikelihoodFunction a, Gen (PrimState m))
  -- ^ Where do we go?
  -> Move m a
  -> m (Link a, LogLikelihoodFunction a, Gen (PrimState m))
-- mhMove = undefined
mhMove ((Link x lX), f, g) (Move p q) = do
  -- 1. Sample new state.
  y <- p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let
    lY = f y
    r  = mhRatio lX lY (q x y) (q y x)
  -- 3. Accept or reject.
  -- XXX: Can this be improved in terms of speed?
  if r >= 0.0
    then return (Link y lY, f, g)
    else do b <- uniform g
            if b < exp (ln r)
            then return (Link y lY, f, g)
            else return (Link x lX, f, g)

mhCycle :: PrimMonad m
  -- ^ Current state and log-likelihood, as well as the log-likelihood function
  -- and a generator.
  => (Link a, LogLikelihoodFunction a, Gen (PrimState m))
  -- ^ Where do we go?
  -> Cycle m a
  -> m (Link a)
-- mhCycle = fold moves in cycle
mhCycle x (Cycle ms) = do
  (l, _, _) <- foldM mhMove x ms
  return l

-- | Run a given number of Metropolis-Hastings cycles.
--
-- TODO: Improve documentation.
mh :: PrimMonad m => Int -> Mcmc m a -> Gen (PrimState m) -> m (Chain a)
-- mh n = repeat mhCycle n times
mh = undefined
