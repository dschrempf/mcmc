{-# LANGUAGE BangPatterns #-}

{- |
Module      :  Mcmc.Metropolis
Description :  Metropolis-Hastings at its best
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 20:11:30 2020.

-}

-- TODO: Move all we can into Mcmc.Mcmc, because it will be reused across
-- algorithms.

module Mcmc.Metropolis
  ( mh
  )
where

import           Prelude                 hiding ( cycle )

import           Control.Monad
import           Control.Monad.IO.Class
import           Control.Monad.Trans.State.Strict
import           Data.Aeson
import           Numeric.Log
import           System.Random.MWC

import           Mcmc.Item
import           Mcmc.Mcmc
import           Mcmc.Move
import           Mcmc.Status
import           Mcmc.Trace

import           Mcmc.Tools.Shuffle

mhRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
mhRatio lX lY qXY qYX = lY * qYX / lX / qXY
{-# INLINE mhRatio #-}

mhMove :: Move a -> Mcmc a ()
mhMove m = do
  let p = mvSample $ mvSimple m
      q = mvLogDensity $ mvSimple m
  s <- get
  let (Item x pX lX) = item s
      pF             = logPriorF s
      lF             = logLikelihoodF s
      a              = acceptance s
      g              = generator s
  -- 1. Sample new state.
  !y <- liftIO $ p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let !pY = pF y
      !lY = lF y
      !r  = mhRatio (pX * lX) (pY * lY) (q x y) (q y x)
  -- 3. Accept or reject.
  if ln r >= 0.0
    then put $ s { item = Item y pY lY, acceptance = prependA m True a }
    else do
      b <- uniform g
      if b < exp (ln r)
        then put $ s { item = Item y pY lY, acceptance = prependA m True a }
        else put $ s { acceptance = prependA m False a }

-- Replicate 'Move's according to their weights and shuffle them.
getNCycles :: Cycle a -> Int -> GenIO -> IO [[Move a]]
getNCycles c = shuffleN mvs
  where !mvs = concat [ replicate (mvWeight m) m | m <- fromCycle c ]

-- Run one iterations; perform all moves in a Cycle.
mhIter :: ToJSON a => [Move a] -> Mcmc a ()
mhIter mvs = do
  mcmcMonitorExec
  mapM_ mhMove mvs
  s <- get
  let i   = item s
      t   = trace s
      n   = iteration s
  put $ s { trace = pushT i t, iteration = succ n }

-- Run N iterations.
mhNIter :: ToJSON a => Int -> Mcmc a ()
mhNIter n = do
  c      <- gets cycle
  g      <- gets generator
  cycles <- liftIO $ getNCycles c n g
  forM_ cycles mhIter

mhBurn :: ToJSON a => Int -> Mcmc a ()
mhBurn n = do
  mhNIter n
  mcmcAutotune n

-- TODO: This is a hack; mhRun should be split into mhBurnIn, and mhRun.
mhBurnIn :: ToJSON a => Int -> Maybe Int -> Mcmc a ()
mhBurnIn b (Just t)
  | t <= 0    = error "mhBurnIn: Auto tuning period smaller equal 0."
  | b > t     = mhBurn t >> mhBurnIn (b - t) (Just t)
  | b <= t    = mhBurn b
  | otherwise = error "mhBurnIn: Please contact maintainer."
mhBurnIn b Nothing = mhNIter b

-- TODO: This is a hack; mhRun should be split into mhBurnIn, and mhRun.
mhRun :: ToJSON a => Maybe Int -> Maybe Int -> Int -> Mcmc a ()
mhRun (Just b) t n
  | b <= 0 = error "mhRun: Number of burn in iterations smaller equal 0."
  | otherwise = do
    liftIO $ putStrLn $ "-- Burn in for " <> show b <> " cycles."
    mcmcMonitorHeader
    mhBurnIn b t
    liftIO $ putStrLn "-- Burn in finished."
    case t of
      Nothing -> return ()
      Just _  -> mcmcSummarizeCycle t
    s <- get
    let a = acceptance s
    put $ s { acceptance = resetA a}
    mhRun Nothing Nothing n
mhRun Nothing _ n = do
  liftIO $ putStrLn $ "-- Run chain for " <> show n <> " iterations."
  mcmcMonitorHeader
  mhNIter n

-- | Run a Markov chain for a given number of Metropolis-Hastings steps.
--
-- Of course, the given status can also be the result of a paused chain.
mh
  :: ToJSON a
  => Status a  -- ^ Initial (or last) 'Status' of the Markov chain.
  -> IO (Status a)
mh = execStateT
  (do
    mcmcInit
    mcmcReport
    mcmcSummarizeCycle Nothing
    -- TODO: This is a hack; mhRun should be split into mhBurnIn, and mhRun.
    s <- get
    let b = burnInIterations s
        t = autoTuningPeriod s
        n = iterations s
    mhRun b t n
    mcmcMonitorExec
    mcmcSummarizeCycle (Just n)
    liftIO $ putStrLn "-- Metropolis-Hastings sampler finished."
    mcmcClose
  )
