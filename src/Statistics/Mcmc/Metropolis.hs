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
  let p = mvSample $ mvSimple m
      q = mvLogDensity $ mvSimple m
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
    then put s { item = Item y lY, acceptance = prependA m True a }
    else do b <- uniform g
            -- Only update the 'Item' after a full cycle.
            if b < exp (ln r)
            then put s { item = Item y lY, acceptance = prependA m True a }
            else put s { acceptance = prependA m False a }

-- Replicate 'Move's according to their weights and shuffle them.
getNCycles :: Cycle a -> Int -> GenIO -> IO [[Move a]]
getNCycles c = shuffleN mvs
  where mvs = concat [ replicate (mvWeight m) m | m <- fromCycle c ]

-- Run a cycle.
mhCycle :: [Move a] -> Mcmc a ()
mhCycle mvs = do
  mapM_ mhMove mvs
  s <- get
  let i = item s
      t = trace s
      n = iteration s
      s' = s { trace = prependT i t, iteration = succ n }
  put s'
  mcmcExecMonitors

-- Run N cycles.
mhNCycles :: Int -> Mcmc a ()
mhNCycles n = do
  c <- gets cycle
  g <- gets generator
  cycles <- liftIO $ getNCycles c n g
  forM_ cycles mhCycle

mhBurn :: Int -> Mcmc a ()
mhBurn n = do
  mhNCycles n
  mcmcAutotune n

-- Burn in with or without auto tuning.
mhBurnIn :: Int -> Maybe Int -> Mcmc a ()
mhBurnIn b (Just t)
  | t <= 0    = error "mhBurnIn: Auto tune period smaller equal 0."
  | b >  t    = mhBurn t >> mhBurnIn (b-t) (Just t)
  | b <= t    = mhBurn b
  | otherwise = error "mhRun: Please contact maintainer."
mhBurnIn b  Nothing = mhNCycles b

mhRun :: Maybe Int -> Maybe Int -> Int -> Mcmc a ()
mhRun (Just b) t n
  | b <= 0    = error "mhBurnIn: Number of burn in cycles smaller equal 0."
  | otherwise = do
      liftIO $ putStrLn $ "-- Burn in for " <> show b <> " cycles."
      mhBurnIn b t
      liftIO $ putStrLn   "-- Burn in finished."
      mhRun Nothing Nothing n
mhRun Nothing _ n = do
  liftIO $ putStrLn $ "-- Start run with " <> show n <> " cycles."
  mhNCycles n

-- | Run a Markov chain for a given number of Metropolis-Hastings steps.
--
-- Of course, the given status can also be the result of a paused chain.
mh :: Show a
  => Maybe Int -- ^ Number of burn in cycles.
  -> Maybe Int -- ^ Auto tune period (only during burn in); deactivate auto
               -- tuning completely with 'Nothing'.
  -> Int       -- ^ Number of Metropolis-Hastings cycles (without auto tuning).
  -> Status a  -- ^ Initial state of Markov chain.
  -> IO (Status a)
mh b t n s =
  execStateT (do liftIO $ putStrLn "-- Start of Metropolis-Hastings sampler."
                 liftIO $ putStr $ summarizeCycle (cycle s)
                 mcmcOpenMonitors
                 liftIO $ putStrLn "-- Initial state of chain."
                 mcmcExecMonitors
                 mhRun b t n
                 mcmcCloseMonitors
                 liftIO $ putStrLn "-- Metropolis-Hastings sampler finished.") s
