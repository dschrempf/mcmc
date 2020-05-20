{- |
Module      :  Statistics.Mcmc.Metropolis
Description :  Metropolis-Hastings at its best
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 20:11:30 2020.

-}

module Statistics.Mcmc.Metropolis
  ( mh
  ) where

import Prelude hiding (cycle)

import Control.Monad
import Control.Monad.Primitive
import Control.Monad.Trans.State.Strict
import Numeric.Log
import System.Random.MWC

import Statistics.Mcmc.Acceptance
import Statistics.Mcmc.Item
import Statistics.Mcmc.Move.Types
import Statistics.Mcmc.Status
import Statistics.Mcmc.Trace

import Statistics.Mcmc.Tools.Shuffle

-- import Debug.Trace

mhRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
mhRatio lX lY qXY qYX = lY * qYX / lX / qXY
{-# INLINE mhRatio #-}

mhMove :: PrimBase m => Move m a -> Mcmc m a ()
mhMove m@(Move _ p q) = do
  s <- get
  let (Item x lX) = item s
      f           = logPosteriorF s
      g           = generator s
      a           = acceptance s
  -- 1. Sample new state.
  y <- liftPrim $ p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let
    lY = f y
    r  = mhRatio lX lY (q x y) (q y x)
  -- 3. Accept or reject.
  -- XXX: Can this be improved in terms of speed?
  -- if traceShow (lX, lY, q x y, q y x, r) $ r >= 1.0
  if ln r >= 0.0
    then put s{item = Item y lY, acceptance = prependA m True a}
    else do b <- uniform g
            -- Only update the 'Item' after a full cycle.
            if b < exp (ln r)
            then put s{item = Item y lY, acceptance = prependA m True a}
            else put s{acceptance = prependA m False a}

-- Replicate 'Move's according to their weights and shuffle them.
shuffleCycle :: PrimMonad m => Cycle m a -> Gen (PrimState m) -> m [Move m a]
shuffleCycle c = shuffle mvs
  where mvs = concat [ replicate w m | (m, w) <- fromCycle c ]

mhCycle :: PrimBase m => Mcmc m a ()
mhCycle = do
  c <- cycle <$> get
  g <- generator <$> get
  mvs <- liftPrim $ shuffleCycle c g
  mapM_ mhMove mvs
  s <- get
  let i = item s
      t = trace s
      s' = s {trace = prependT i t}
  put s'

-- Run a given number of Metropolis-Hastings cycles.
mhRun :: PrimBase m => Int -> Mcmc m a ()
mhRun n = replicateM_ n mhCycle

-- | Run a Markov chain for a given number of Metropolis-Hastings steps.
--
-- Of course, the given status can also be the result of a paused chain.
mh :: (Show a, PrimBase m)
  => Int -- ^ Number of Metropolis-Hastings steps.
  -> Status m a -- ^ Initial state of Markov chain.
  -> m (Status m a)
mh n = execStateT (mhRun n)
