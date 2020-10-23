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
--
-- Metropolis-Hastings algorithm.
module Mcmc.Metropolis
  ( mh,
    mhContinue,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State
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

-- The Metropolis-Hastings ratio.
--
-- 'Infinity' if fX is zero. In this case, the proposal is always accepted.
--
-- 'NaN' if (fY or q) and fX are zero. In this case, the proposal is always
-- rejected.

-- There is a discrepancy between authors saying that one should (a) always
-- accept the new state when the current posterior is zero (Chapter 4 of the
-- Handbook of Markov Chain Monte Carlo), or (b) almost surely reject the
-- proposal when either fY or q are zero (Chapter 1). Since I trust the author
-- of Chapter 1 (Charles Geyer) I choose to follow option (b).
mhRatio :: Log Double -> Log Double -> Log Double -> Log Double
-- q = qYX / qXY
mhRatio fX fY q = fY * q / fX
{-# INLINE mhRatio #-}

mhPropose :: Proposal a -> Mcmc a ()
mhPropose m = do
  let p = pSimple m
  s <- get
  let (Item x pX lX) = item s
      pF = priorF s
      lF = likelihoodF s
      a = acceptance s
      g = generator s
  -- 1. Sample new state.
  (!y, !q) <- liftIO $ p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let !pY = pF y
      !lY = lF y
      !r = mhRatio (pX * lX) (pY * lY) q
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
  mcmcClean
  mcmcMonitorExec

-- Run N iterations.
mhNIter :: ToJSON a => Int -> Mcmc a ()
mhNIter n = do
  mcmcDebugS $ "Run " <> show n <> " iterations."
  c <- gets cycle
  g <- gets generator
  cycles <- liftIO $ getNIterations c n g
  forM_ cycles mhIter

-- Burn in and auto tune.
mhBurnInN :: ToJSON a => Int -> Maybe Int -> Mcmc a ()
mhBurnInN b (Just t)
  | t <= 0 = error "mhBurnInN: Auto tuning period smaller equal 0."
  | b > t = do
    mcmcResetA
    mhNIter t
    mcmcSummarizeCycle >>= mcmcDebugB
    mcmcAutotune
    mhBurnInN (b - t) (Just t)
  | otherwise = do
    mcmcResetA
    mhNIter b
    mcmcSummarizeCycle >>= mcmcInfoB
    mcmcInfoS $ "Acceptance ratios calculated over the last " <> show b <> " iterations."
mhBurnInN b Nothing = mhNIter b

-- Initialize burn in for given number of iterations.
mhBurnIn :: ToJSON a => Int -> Maybe Int -> Mcmc a ()
mhBurnIn b t
  | b < 0 = error "mhBurnIn: Negative number of burn in iterations."
  | b == 0 = return ()
  | otherwise = do
    mcmcInfoS $ "Burn in for " <> show b <> " cycles."
    mcmcDebugS $ "Auto tuning period is " <> show t <> "."
    mhBurnInN b t
    mcmcInfoB "Burn in finished."

-- Run for given number of iterations.
mhRun :: ToJSON a => Int -> Mcmc a ()
mhRun n = do
  mcmcResetA
  mcmcInfoS $ "Run chain for " <> show n <> " iterations."
  -- let (m, r) = n `quotRem` 100
  -- -- Print header to standard output every 100 iterations.
  -- replicateM_ m $ do
  --   mcmcMonitorStdOutHeader
  --   mhNIter 100
  -- when (r > 0) $ do
  --   mcmcMonitorStdOutHeader
  --   mhNIter r
  mhNIter n

mhT :: ToJSON a => Mcmc a ()
mhT = do
  mcmcInfoB "Metropolis-Hastings sampler."
  mcmcSummarizeCycle >>= mcmcInfoB
  mcmcReport
  s <- get
  let b = fromMaybe 0 (burnInIterations s)
  mhBurnIn b (autoTuningPeriod s)
  mhRun $ iterations s

mhContinueT :: ToJSON a => Int -> Mcmc a ()
mhContinueT dn = do
  mcmcInfoB "Continuation of Metropolis-Hastings sampler."
  mcmcInfoS $ "Run chain for " <> show dn <> " additional iterations."
  mcmcSummarizeCycle >>= mcmcInfoB
  mhRun dn

-- | Continue a Markov chain for a given number of Metropolis-Hastings steps.
--
-- At the moment, when an MCMC run is continued, the old @.mcmc@ file is
-- deleted. This behavior may change in the future.
--
-- This means that an interrupted continuation also breaks previous runs. This
-- step is necessary because, otherwise, incomplete monitor files are left on
-- disk, if a continuation is canceled. Subsequent continuations would append to
-- the incomplete monitor files and produce garbage.
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
  where
    n' = iterations s + dn
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
