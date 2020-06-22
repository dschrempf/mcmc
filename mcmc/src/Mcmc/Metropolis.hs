{-# LANGUAGE BangPatterns #-}

-- |
-- Module      :  Mcmc.Metropolis
-- Description :  Metropolis-Hastings at its best
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 20:11:30 2020.
module Mcmc.Metropolis
  ( mh,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict
import Data.Aeson
import Mcmc.Item
import Mcmc.Mcmc
import Mcmc.Move
import Mcmc.Status
import Mcmc.Tools.Shuffle
import Mcmc.Trace
import Numeric.Log
import System.Random.MWC
import Prelude hiding (cycle)

-- For non-symmetric moves.
mhRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
mhRatio lX lY qXY qYX = lY * qYX / lX / qXY
{-# INLINE mhRatio #-}

-- For symmetric moves.
mhRatioSymmetric :: Log Double -> Log Double -> Log Double
mhRatioSymmetric lX lY = lY / lX
{-# INLINE mhRatioSymmetric #-}

mhMove :: Move a -> Mcmc a ()
mhMove m = do
  let p = mvSample $ mvSimple m
      mq = mvLogDensity $ mvSimple m
  s <- get
  let (Item x pX lX) = item s
      pF = logPriorF s
      lF = logLikelihoodF s
      a = acceptance s
      g = generator s
  -- 1. Sample new state.
  !y <- liftIO $ p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let !pY = pF y
      !lY = lF y
      !r = case mq of
        Nothing -> mhRatioSymmetric (pX * lX) (pY * lY)
        Just q -> mhRatio (pX * lX) (pY * lY) (q x y) (q y x)
  -- 3. Accept or reject.
  if ln r >= 0.0
    then put $ s {item = Item y pY lY, acceptance = pushA m True a}
    else do
      b <- uniform g
      if b < exp (ln r)
        then put $ s {item = Item y pY lY, acceptance = pushA m True a}
        else put $ s {acceptance = pushA m False a}

-- Replicate 'Move's according to their weights and shuffle them.
getNCycles :: Cycle a -> Int -> GenIO -> IO [[Move a]]
getNCycles c = shuffleN mvs
  where
    !mvs = concat [replicate (mvWeight m) m | m <- fromCycle c]

-- Run one iterations; perform all moves in a Cycle.
mhIter :: ToJSON a => [Move a] -> Mcmc a ()
mhIter mvs = do
  mcmcMonitorExec
  mapM_ mhMove mvs
  s <- get
  let i = item s
      t = trace s
      n = iteration s
  put $ s {trace = pushT i t, iteration = succ n}

-- Run N iterations.
mhNIter :: ToJSON a => Int -> Mcmc a ()
mhNIter n = do
  c <- gets cycle
  g <- gets generator
  cycles <- liftIO $ getNCycles c n g
  forM_ cycles mhIter

-- Burn in and auto tune.
mhBurnInN :: ToJSON a => Int -> Maybe Int -> Mcmc a ()
mhBurnInN b (Just t)
  | t <= 0 = error "mhBurnInN: Auto tuning period smaller equal 0."
  | b > t = mhNIter t >> mcmcAutotune t >> mhBurnInN (b - t) (Just t)
  | otherwise = mhNIter b >> mcmcAutotune b
mhBurnInN b Nothing = mhNIter b

-- Initialize burn in for given number of iterations.
mhBurnIn :: ToJSON a => Int -> Maybe Int -> Mcmc a ()
mhBurnIn b t
  | b <= 0 = error "mhBurnIn: Number of burn in iterations smaller equal 0."
  | otherwise = do
    liftIO $ putStrLn $ "-- Burn in for " <> show b <> " cycles."
    mcmcMonitorHeader
    mhBurnInN b t
    liftIO $ putStrLn "-- Burn in finished."
    case t of
      Nothing -> return ()
      Just _ -> mcmcSummarizeCycle t
    s <- get
    let a = acceptance s
    put $ s {acceptance = resetA a}

-- Run for given number of iterations.
mhRun :: ToJSON a => Int -> Mcmc a ()
mhRun n = do
  liftIO $ putStrLn $ "-- Run chain for " <> show n <> " iterations."
  mcmcMonitorHeader
  mhNIter n

mhT :: ToJSON a => Mcmc a ()
mhT = do
  s <- get
  if iteration s == 0
    then liftIO $ putStrLn "-- Continue Metropolis-Hastings sampler."
    else liftIO $ putStrLn "-- Start of Metropolis-Hastings sampler."
  mcmcInit
  mcmcReport
  mcmcSummarizeCycle Nothing
  case burnInIterations s of
    Nothing -> return ()
    Just b -> mhBurnIn b (autoTuningPeriod s)
  mhRun (iterations s)
  mcmcMonitorExec
  mcmcSummarizeCycle (Just $ iterations s)
  liftIO $ putStrLn "-- Metropolis-Hastings sampler finished."
  mcmcClose

-- | Run a Markov chain for a given number of Metropolis-Hastings steps.
--
-- Of course, the given status can also be the result of a paused chain.
mh ::
  ToJSON a =>
  -- | Initial (or last) status of the Markov chain.
  Status a ->
  IO (Status a)
mh = execStateT mhT
