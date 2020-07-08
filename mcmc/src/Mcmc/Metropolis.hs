{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

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
import Mcmc.Proposal
import Mcmc.Status
import Mcmc.Trace
import Numeric.Log
import System.Random.MWC
import Prelude hiding (cycle)

-- For non-symmetric proposals.
mhRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
mhRatio fX fY qXY qYX = fY * qYX / fX / qXY
{-# INLINE mhRatio #-}

-- For symmetric proposals.
mhRatioSymmetric :: Log Double -> Log Double -> Log Double
mhRatioSymmetric fX fY = fY / fX
{-# INLINE mhRatioSymmetric #-}

mhPropose :: Proposal a -> Mcmc a ()
mhPropose m = do
  let p = pSample $ pSimple m
      mq = pKernel $ pSimple m
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

-- Run one iterations; perform all proposals in a Cycle.
mhIter :: ToJSON a => [Proposal a] -> Mcmc a ()
mhIter ps = do
  mapM_ mhPropose ps
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
  | b > t = do
    mcmcResetA
    mhNIter t
    mcmcDebug mcmcSummarizeCycle
    mcmcAutotune
    mhBurnInN (b - t) (Just t)
  | otherwise = do
    mcmcResetA
    mhNIter b
    mcmcInfo mcmcSummarizeCycle
    mcmcInfo $ mcmcOutS $ "Acceptance ratios calculated over the last " <> show b <> " iterations."
mhBurnInN b Nothing = mhNIter b

-- Initialize burn in for given number of iterations.
mhBurnIn :: ToJSON a => Int -> Maybe Int -> Mcmc a ()
mhBurnIn b t
  | b < 0 = error "mhBurnIn: Negative number of burn in iterations."
  | b == 0 = return ()
  | otherwise = do
    mcmcInfo $ mcmcOutS $ "Burn in for " <> show b <> " cycles."
    mcmcMonitorStdOutHeader
    mhBurnInN b t
    mcmcInfo $ mcmcOut "Burn in finished."

-- Run for given number of iterations.
mhRun :: ToJSON a => Int -> Mcmc a ()
mhRun n = do
  mcmcInfo $ mcmcOutS $ "Run chain for " <> show n <> " iterations."
  mcmcMonitorStdOutHeader
  mhNIter n

mhT :: ToJSON a => Mcmc a ()
mhT = do
  mcmcInfo $ mcmcOut "-- Metropolis-Hastings sampler."
  mcmcReport
  mcmcInfo mcmcSummarizeCycle
  s <- get
  let b = fromMaybe 0 (burnInIterations s)
  mhBurnIn b (autoTuningPeriod s)
  let n = iterations s
  mhRun n

mhContinueT :: ToJSON a => Int -> Mcmc a ()
mhContinueT dn = do
  mcmcInfo $ mcmcOut "-- Continuation of Metropolis-Hastings sampler."
  mcmcInfo $ mcmcOutS $ "-- Run chain for " <> show dn <> " additional iterations."
  mcmcInfo mcmcSummarizeCycle
  mhRun dn

-- | Continue a Markov chain for a given number of Metropolis-Hastings steps.
mhContinue ::
  ToJSON a =>
  -- | Additional number of Metropolis-Hastings steps.
  Int ->
  -- | Loaded status of the Markov chain.
  Status a ->
  IO (Status a)
mhContinue dn s
  | dn <= 0 = error "mhContinue: The number of iterations is zero or negative."
  | otherwise = mcmcRun (mhContinueT dn) s'
    where n' = iterations s + dn
          s' = s {iterations = n'}

-- | Run a Markov chain for a given number of Metropolis-Hastings steps.
mh ::
  ToJSON a =>
  -- | Initial (or last) status of the Markov chain.
  Status a ->
  IO (Status a)
mh s =
  if iteration s == 0
    then mcmcRun mhT s
    else do
      putStrLn "To continue a Markov chain run, please use 'mhContinue'."
      error $ "mh: Current iteration " ++ show (iteration s) ++ " is non-zero."
