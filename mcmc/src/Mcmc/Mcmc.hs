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
--
-- Run an MCMC algorithm.
module Mcmc.Mcmc
  ( mcmc,
    mcmcContinue,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.RWS.CPS
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Time.Clock
import Data.Time.Format
import Mcmc.Algorithm
import Mcmc.Environment
import Mcmc.Monitor.Time
import Mcmc.Settings
import System.IO
import Prelude hiding (cycle)

-- The MCMC algorithm has read access to an environment and uses an algorithm
-- transforming the state @a@.
type MCMC a = RWST Environment () a IO

msgPrepare :: Char -> BL.ByteString -> BL.ByteString
msgPrepare c t = BL.cons c $ ": " <> t

-- Write to standard output and log file.
mcmcOutB :: BL.ByteString -> MCMC a ()
mcmcOutB msg = do
  h <- reader logHandle
  liftIO $ BL.putStrLn msg >> BL.hPutStrLn h msg

-- -- Perform warning action.
-- mcmcWarnA :: MCMC a () -> MCMC a ()
-- mcmcWarnA a = reader (verbosity . settings) >>= \v -> when (v >= Warn) a

-- -- Print warning message.
-- mcmcWarnB :: BL.ByteString -> MCMC a ()
-- mcmcWarnB = mcmcWarnA . mcmcOutB . msgPrepare 'W'

-- -- Print warning message.
-- mcmcWarnS :: String -> MCMC a ()
-- mcmcWarnS = mcmcWarnB . BL.pack

-- Perform info action.
mcmcInfoA :: MCMC a () -> MCMC a ()
mcmcInfoA a = reader (sVerbosity . settings) >>= \v -> when (v >= Info) a

-- Print info message.
mcmcInfoB :: BL.ByteString -> MCMC a ()
mcmcInfoB = mcmcInfoA . mcmcOutB . msgPrepare 'I'

-- Print info message.
mcmcInfoS :: String -> MCMC a ()
mcmcInfoS = mcmcInfoB . BL.pack

-- Perform debug action.
mcmcDebugA :: MCMC a () -> MCMC a ()
mcmcDebugA a = reader (sVerbosity . settings) >>= \v -> when (v == Debug) a

-- Print debug message.
mcmcDebugB :: BL.ByteString -> MCMC a ()
mcmcDebugB = mcmcDebugA . mcmcOutB . msgPrepare 'D'

-- Print debug message.
mcmcDebugS :: String -> MCMC a ()
mcmcDebugS = mcmcDebugB . BL.pack

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

mcmcReportTime :: MCMC a ()
mcmcReportTime = do
  ti <- reader startingTime
  mcmcInfoS $ "Starting time: " <> fTime ti

mcmcExecute :: Algorithm t a => MCMC (t a) ()
mcmcExecute = do
  s <- reader settings
  case sExecutionMode s of
    Fail -> mcmcNewRun
    Overwrite -> mcmcNewRun
    Continue -> mcmcContinueRun

mcmcNewRun :: Algorithm t a => MCMC (t a) ()
mcmcNewRun = do
  s <- reader settings
  mcmcInfoB "Start new MCMC sampler."
  mcmcInfoB "Initial state."
  mcmcExecuteMonitors
  get >>= mcmcInfoB . aSummarizeCycle
  mcmcBurnIn
  mcmcResetAcceptance
  let i = sIterations s
  mcmcInfoS $ "Run chain for " ++ show i ++ " iterations."
  mcmcIterate i

mcmcContinueRun :: Algorithm t a => MCMC (t a) ()
mcmcContinueRun = do
  s <- reader settings
  let iTotal = sIterations s + burnInIterations (sBurnIn s)
  mcmcInfoB "Continuation of MCMC sampler."
  a <- get
  let iCurrent = aIteration a
  mcmcInfoS $ "Current iteration: " ++ show iCurrent ++ "."
  mcmcInfoS $ "Total iterations: " ++ show iTotal ++ "."
  let di = iTotal - iCurrent
  get >>= mcmcInfoB . aSummarizeCycle
  mcmcInfoS $ "Run chain for " ++ show di ++ " iterations."
  mcmcIterate di

mcmcBurnIn :: Algorithm t a => MCMC (t a) ()
mcmcBurnIn = do
  s <- reader settings
  case sBurnIn s of
    NoBurnIn ->
      mcmcInfoS "No burn in."
    BurnInNoAutoTuning n -> do
      mcmcInfoS $ "Burn in for " <> show n <> " iterations."
      mcmcInfoS "Auto tuning is disabled."
      mcmcIterate n
      get >>= mcmcInfoB . aSummarizeCycle
      mcmcInfoB "Burn in finished."
    BurnInWithAutoTuning n t -> do
      mcmcInfoS $ "Burn in for " ++ show n ++ " iterations."
      mcmcInfoS $ "Auto tuning is enabled with a period of " ++ show t ++ "."
      mcmcBurnInWithAutoTuning n t
      mcmcInfoB "Burn in finished."

mcmcBurnInWithAutoTuning :: Algorithm t a => Int -> Int -> MCMC (t a) ()
mcmcBurnInWithAutoTuning b t
  | b > t = do
    mcmcResetAcceptance
    mcmcIterate t
    get >>= mcmcDebugB . aSummarizeCycle
    mcmcAutotune
    mcmcBurnInWithAutoTuning (b - t) t
  | otherwise = do
    mcmcResetAcceptance
    mcmcIterate b
    get >>= mcmcInfoB . aSummarizeCycle
    mcmcInfoS $ "Acceptance ratios calculated over the last " <> show b <> " iterations."

mcmcIterate :: Algorithm t a => Int -> MCMC (t a) ()
mcmcIterate n
  | n < 0 = error "mcmcIterate: Number of iterations is negative."
  | n == 0 = return ()
  | otherwise = do
    -- TODO: Splitmix. Remove IO monad as soon as possible.
    get >>= liftIO . aIterate >>= put
    mcmcExecuteMonitors
    mcmcIterate (n -1)

-- Execute the monitors of the chain.
mcmcExecuteMonitors :: Algorithm t a => MCMC (t a) ()
mcmcExecuteMonitors = do
  e <- ask
  a <- get
  mStdLog <- liftIO (aExecuteMonitors e a)
  case mStdLog of
    Nothing -> return ()
    Just x -> mcmcOutB x

-- Auto tune the proposals.
mcmcAutotune :: Algorithm t a => MCMC (t a) ()
mcmcAutotune = do
  mcmcDebugB "Auto tune."
  modify aAutoTune

-- Reset acceptance counts.
mcmcResetAcceptance :: Algorithm t a => MCMC (t a) ()
mcmcResetAcceptance = do
  mcmcDebugB "Reset acceptance ratios."
  modify aResetAcceptance

-- Report and finish up.
mcmcClose :: Algorithm t a => MCMC (t a) ()
mcmcClose = do
  a <- get
  mcmcInfoB $ aSummarizeCycle a
  mcmcInfoS $ aName a ++ " algorithm finished."
  mcmcSave
  ti <- reader startingTime
  te <- liftIO getCurrentTime
  let dt = te `diffUTCTime` ti
  mcmcInfoB $ "Wall clock run time: " <> renderDuration dt <> "."
  mcmcInfoS $ "End time: " <> fTime te
  h <- reader logHandle
  liftIO $ hClose h

-- Save the MCMC run.
mcmcSave :: Algorithm t a => MCMC (t a) ()
mcmcSave = do
  ss <- reader settings
  a <- get
  case sSaveMode ss of
    NoSave -> mcmcInfoB "Do not save the Markov chain."
    SaveWithTrace n -> do
      mcmcInfoB "Save settings."
      liftIO $ settingsSave ss
      let nm = sAnalysisName ss
      mcmcInfoS $
        "Save compressed Markov chain with trace of length " ++ show n ++ "."
      mcmcInfoB "For long traces, or complex objects, this may take a while."
      liftIO $ aSave n nm a
      mcmcInfoB "Markov chain saved."

-- Initialize the run, execute the run, and close the run.
mcmcRun :: Algorithm t a => MCMC (t a) ()
mcmcRun = do
  mcmcDebugB "The settings are:"
  reader settings >>= mcmcDebugS . show

  -- Initialize.
  a <- get
  mcmcInfoS $ aName a ++ " algorithm."
  e <- ask
  get >>= liftIO . aOpenMonitors e >>= put
  mcmcReportTime

  -- Execute.
  mcmcExecute

  -- Close.
  mcmcClose
  get >>= liftIO . aCloseMonitors >>= put

-- | Run an MCMC algorithm with given settings.
mcmc :: Algorithm t a => Settings -> t a -> IO (t a)
mcmc s a = do
  settingsCheck s $ aIteration a
  e <- initializeEnvironment s
  fst <$> execRWST mcmcRun e a

-- | Continue an MCMC algorithm for the given number of iterations.
mcmcContinue :: Algorithm t a => Int -> Settings -> t a -> IO (t a)
mcmcContinue dn s = mcmc s'
  where
    n = sIterations s
    s' = s {sIterations = n + dn, sExecutionMode = Continue}
