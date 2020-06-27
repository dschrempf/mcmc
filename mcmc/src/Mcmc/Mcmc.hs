{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Mcmc
-- Description :  Mcmc helpers
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 10:19:45 2020.
module Mcmc.Mcmc
  ( Mcmc,
    mcmcInfo,
    mcmcAutotune,
    mcmcResetA,
    mcmcSummarizeCycle,
    mcmcInit,
    mcmcReport,
    mcmcMonitorHeader,
    mcmcMonitorExec,
    mcmcClose,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict
import Data.Aeson
import Data.Maybe
import qualified Data.Text.Lazy.IO as T
import Data.Time.Clock
import Data.Time.Format
import Mcmc.Monitor
import Mcmc.Monitor.Time
import Mcmc.Move
import Mcmc.Save
import Mcmc.Status
import Mcmc.Verbosity
import Prelude hiding (cycle)

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO

-- -- Print warning.
-- mWarn :: IO () -> Mcmc a ()
-- mWarn m = gets verbosity >>= \v -> liftIO (warn v m)

-- | Perform IO action only if 'Verbosity' is 'Info' or larger.
mcmcInfo :: IO () -> Mcmc a ()
mcmcInfo m = gets verbosity >>= \v -> liftIO (info v m)

-- -- Print debugging message.
-- mDebug :: IO () -> Mcmc a ()
-- mDebug m = gets verbosity >>= \v -> liftIO (debug v m)

-- | Auto tune the 'Move's in the 'Cycle' of the chain. Reset acceptance counts.
-- See 'autotune'.
mcmcAutotune :: Mcmc a ()
mcmcAutotune = do
  s <- get
  let a = acceptance s
      c = cycle s
      c' = autotuneCycle a c
  put $ s {cycle = c'}

-- | Reset acceptance counts.
mcmcResetA :: Mcmc a ()
mcmcResetA = do
  s <- get
  let a = acceptance s
  put $ s {acceptance = resetA a}

-- TODO: This function should not have a Bool parameter. Rather, the
-- summarizeCycle function should be aware of an empty acceptance type.

-- | Print short summary of 'Move's in 'Cycle'. If 'True', also print acceptance
-- ratios. See 'summarizeCycle'.
mcmcSummarizeCycle :: Bool -> Mcmc a ()
mcmcSummarizeCycle aPred = do
  a <- gets acceptance
  let ma = if aPred then Just a else Nothing
  c <- gets cycle
  mcmcInfo $ T.putStr $ summarizeCycle ma c

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

-- | Set the total number of iterations, the current time and open the
-- 'Monitor's of the chain. See 'mOpen'.
mcmcInit :: Mcmc a ()
mcmcInit = do
  t <- liftIO getCurrentTime
  mcmcInfo $ putStrLn $ "-- Start time: " <> fTime t
  s <- get
  let m = monitor s
      n = iteration s
      nm = name s
  m' <-
    if n == 0
      then liftIO $ mOpen nm m
      else liftIO $ mAppend nm m
  put $ s {monitor = m', start = Just (n, t)}

-- | Report what is going to be done.
mcmcReport :: ToJSON a => Mcmc a ()
mcmcReport = do
  s <- get
  let b = burnInIterations s
      t = autoTuningPeriod s
      n = iterations s
  case b of
    Just b' -> mcmcInfo $ putStrLn $ "-- Burn in for " <> show b' <> " iterations."
    Nothing -> return ()
  case t of
    Just t' ->
      mcmcInfo $ putStrLn $
        "-- Auto tune every "
          <> show t'
          <> " iterations (during burn in only)."
    Nothing -> return ()
  mcmcInfo $ putStrLn $ "-- Run chain for " <> show n <> " iterations."
  mcmcInfo $ putStrLn "-- Initial state."
  mcmcMonitorHeader
  mcmcMonitorExec

-- | Print header line of 'Monitor' (only standard output).
mcmcMonitorHeader :: Mcmc a ()
mcmcMonitorHeader = do
  m <- gets monitor
  mcmcInfo $ mHeader m

-- Save the status of an MCMC run. See 'saveStatus'.
mcmcSave :: ToJSON a => Mcmc a ()
mcmcSave = do
  s <- get
  if save s
    then do
      mcmcInfo $ putStrLn "-- Save Markov chain. For long chains, this may take a while."
      liftIO $ saveStatus (name s <> ".mcmc") s
      mcmcInfo $ putStrLn "-- Done saving Markov chain."
    else mcmcInfo $ putStrLn "-- Do not save the Markov chain."

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: ToJSON a => Mcmc a ()
mcmcMonitorExec = do
  s <- get
  let i = iteration s
      j = iterations s + fromMaybe 0 (burnInIterations s)
      m = monitor s
      st = fromMaybe (error "mcmcMonitorExec: Starting state and time not set.") (start s)
      tr = trace s
      vb = verbosity s
  liftIO $ mExec vb i st tr j m

-- | Close the 'Monitor's of the chain. See 'mClose'.
mcmcClose :: ToJSON a => Mcmc a ()
mcmcClose = do
  s <- get
  mcmcSummarizeCycle True
  mcmcInfo $ putStrLn "-- Metropolis-Hastings sampler finished."
  let m = monitor s
  m' <- liftIO $ mClose m
  put $ s {monitor = m'}
  mcmcSave
  t <- liftIO getCurrentTime
  let rt = case start s of
        Nothing -> error "mcmcClose: Start time not set."
        Just (_, st) -> t `diffUTCTime` st
  mcmcInfo $ T.putStrLn $ "-- Wall clock run time: " <> renderDuration rt <> "."
  mcmcInfo $ putStrLn $ "-- End time: " <> fTime t
