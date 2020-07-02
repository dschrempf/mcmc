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
    mhContinue,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict
import Data.Aeson
import Data.Maybe
import Mcmc.Item
import Mcmc.Mcmc
import Mcmc.Move
import Mcmc.Status
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
      mq = mvDensity $ mvSimple m
  s <- get
  let (Item x pX lX) = item s
      pF = priorF s
      lF = likelihoodF s
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

-- TODO: Splitmix. Split the generator here. See SaveSpec -> mhContinue.

-- Run one iterations; perform all moves in a Cycle.
mhIter :: ToJSON a => [Move a] -> Mcmc a ()
mhIter mvs = do
  mapM_ mhMove mvs
  s <- get
  let i = item s
      t = trace s
      n = iteration s
  put $ s {trace = pushT i t, iteration = succ n}
  mcmcMonitorExec

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
  | b > t = mhNIter t >> mcmcAutotune >> mcmcDebugA mcmcSummarizeCycle >> mcmcResetA >> mhBurnInN (b - t) (Just t)
  | otherwise = do
      mhNIter b
      mcmcAutotune
      -- TODO: This is not ideal, yet.
      when (b <= 100) (mcmcWarn $ "WARNING: Last auto tuning period spans " <> show b <> " iterations only.")
      mcmcAutotune
      mcmcInfo "Acceptance ratios calculated before the last auto tune."
      mcmcSummarizeCycle
      mcmcResetA
mhBurnInN b Nothing = mhNIter b

-- Initialize burn in for given number of iterations.
mhBurnIn :: ToJSON a => Int -> Maybe Int -> Mcmc a ()
mhBurnIn b t
  | b < 0 = error "mhBurnIn: Negative number of burn in iterations."
  | b == 0 = return ()
  | otherwise = do
    mcmcInfo $ "Burn in for " <> show b <> " cycles."
    mcmcMonitorHeader
    mhBurnInN b t
    mcmcInfo "Burn in finished."

-- Run for given number of iterations.
mhRun :: ToJSON a => Int -> Mcmc a ()
mhRun n = do
  mcmcInfo $ "Run chain for " <> show n <> " iterations."
  mcmcMonitorHeader
  mhNIter n

mhT :: ToJSON a => Mcmc a ()
mhT = do
  s <- get
  mcmcInfo "Start of Metropolis-Hastings sampler."
  mcmcInit
  mcmcReport
  mcmcSummarizeCycle
  let b = fromMaybe 0 (burnInIterations s)
  mhBurnIn b (autoTuningPeriod s)
  let n = iterations s
  mhRun n
  mcmcClose

mhContinueT :: ToJSON a => Int -> Mcmc a ()
mhContinueT dn = do
  mcmcInfo "Continue Metropolis-Hastings sampler."
  mcmcInfo $ "Run chain for " <> show dn <> " additional iterations."
  s <- get
  let n = iterations s
  put s {iterations = n + dn}
  mcmcInit
  mcmcSummarizeCycle
  mhRun dn
  mcmcClose

-- | Continue a Markov chain for a given number of Metropolis-Hastings steps.
mhContinue ::
  ToJSON a =>
  -- | Additional number of Metropolis-Hastings steps.
  Int ->
  -- | Loaded status of the Markov chain.
  Status a ->
  IO (Status a)
mhContinue dn
  | dn <= 0 = error "mhContinue: The number of iterations is zero or negative."
  | otherwise = execStateT $ mhContinueT dn

-- | Run a Markov chain for a given number of Metropolis-Hastings steps.
mh ::
  ToJSON a =>
  -- | Initial (or last) status of the Markov chain.
  Status a ->
  IO (Status a)
mh s = do
  let m = iteration s
  if m == 0
    then execStateT mhT s
    else do
      putStrLn "To continue a Markov chain run, please use 'mhContinue'."
      error $ "mh: Current iteration " ++ show m ++ " is non-zero."
