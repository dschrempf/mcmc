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
import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict
import Numeric.Log
import System.Random.MWC

import Statistics.Mcmc.Acceptance
import Statistics.Mcmc.Item
import Statistics.Mcmc.Move
import Statistics.Mcmc.Status
import Statistics.Mcmc.Trace

import Statistics.Mcmc.Tools.Shuffle

mhRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
mhRatio lX lY qXY qYX = lY * qYX / lX / qXY
{-# INLINE mhRatio #-}

mhMove :: Move a -> Mcmc a ()
mhMove m = do
  let p = mvSample m
      q = mvLogDensity m
  s <- get
  let (Item x lX) = item s
      f           = logPosteriorF s
      g           = generator s
      a           = acceptance s
  -- 1. Sample new state.
  y <- liftIO $ p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let
    lY = f y
    r  = mhRatio lX lY (q x y) (q y x)
  -- 3. Accept or reject.
  if ln r >= 0.0
    then put s{item = Item y lY, acceptance = prependA m True a}
    else do b <- uniform g
            -- Only update the 'Item' after a full cycle.
            if b < exp (ln r)
            then put s{item = Item y lY, acceptance = prependA m True a}
            else put s{acceptance = prependA m False a}

-- Replicate 'Move's according to their weights and shuffle them.
getCycles :: Cycle a -> Int -> GenIO -> IO [[Move a]]
getCycles c = shuffleN mvs
  where mvs = concat [ replicate w m | (m, w) <- fromCycle c ]

mhCycle :: [Move a] -> Mcmc a ()
mhCycle mvs = do
  mapM_ mhMove mvs
  s <- get
  let i = item s
      t = trace s
      n = iteration s
      s' = s {trace = prependT i t, iteration = succ n}
  put s'

-- TODO: Burn in.

-- TODO: Tune.

-- Run a given number of Metropolis-Hastings cycles.
mhRun :: Int -> Mcmc a ()
mhRun n = do
  c <- cycle <$> get
  g <- generator <$> get
  cycles <- liftIO $ getCycles c n g
  forM_ cycles mhCycle

-- | Run a Markov chain for a given number of Metropolis-Hastings steps.
--
-- Of course, the given status can also be the result of a paused chain.
mh :: Show a
  => Int -- ^ Number of Metropolis-Hastings steps.
  -> Status a -- ^ Initial state of Markov chain.
  -> IO (Status a)
mh n = execStateT (mhRun n)
