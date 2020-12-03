{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Mcmc
-- Description :  Framework for running Markov chain Monte Carlo samplers
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 10:19:45 2020.
--
-- This module provides the general framework for running MCMC samplers. By
-- design choice this module is agnostic about the details of the used
-- 'Algorithm'.
module Mcmc.Mcmc
  ( mcmc,
    mcmcContinue,
  )
where

import Control.Monad
import Control.Monad.IO.Class
-- import Control.Monad.Trans.RWS.CPS
import Control.Monad.Trans.Reader
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Maybe
import Data.Time.Clock
import Mcmc.Algorithm
import Mcmc.Environment
import Mcmc.Monitor.Time
import Mcmc.Settings
import System.IO
import Text.Show.Pretty
import Prelude hiding (cycle)

-- The MCMC algorithm has read access to an environment and uses an algorithm
-- transforming the state @a@.
type MCMC a = ReaderT Environment IO a

msgPrepare :: BL.ByteString -> BL.ByteString -> BL.ByteString
msgPrepare pref msg = BL.intercalate "\n" $ map (BL.append pref) $ BL.lines msg

-- Write to standard output and log file.
mcmcOutB :: BL.ByteString -> BL.ByteString -> MCMC ()
mcmcOutB pref msg = do
  h <- fromMaybe (error "mcmcOut: Log handle is missing.") <$> reader logHandle
  liftIO $ BL.putStrLn msg' >> BL.hPutStrLn h msg'
  where
    msg' = msgPrepare pref msg

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
mcmcInfoA :: MCMC () -> MCMC ()
mcmcInfoA a = reader (sVerbosity . settings) >>= \v -> when (v >= Info) a

-- Print info message.
mcmcInfoB :: BL.ByteString -> MCMC ()
mcmcInfoB = mcmcInfoA . mcmcOutB "I: "

-- Print info message.
mcmcInfoS :: String -> MCMC ()
mcmcInfoS = mcmcInfoB . BL.pack

-- Perform debug action.
mcmcDebugA :: MCMC () -> MCMC ()
mcmcDebugA a = reader (sVerbosity . settings) >>= \v -> when (v == Debug) a

-- Print debug message.
mcmcDebugB :: BL.ByteString -> MCMC ()
mcmcDebugB = mcmcDebugA . mcmcOutB "D: "

-- Print debug message.
mcmcDebugS :: String -> MCMC ()
mcmcDebugS = mcmcDebugB . BL.pack

mcmcReportTime :: MCMC ()
mcmcReportTime = do
  mcmcDebugB "Report time."
  ti <- reader startingTime
  mcmcInfoS $ "Starting time of MCMC sampler: " <> renderTime ti

mcmcExecute :: Algorithm a => a -> MCMC a
mcmcExecute a = do
  mcmcDebugB "Executing MCMC run."
  s <- reader settings
  a' <- case sExecutionMode s of
    Fail -> mcmcNewRun a
    Overwrite -> mcmcNewRun a
    Continue -> mcmcContinueRun a
  mcmcDebugB "Executed MCMC run."
  return a'

-- Reset acceptance counts.
mcmcResetAcceptance :: Algorithm a => a -> MCMC a
mcmcResetAcceptance a = do
  mcmcDebugB "Reset acceptance rates."
  return $ aResetAcceptance a

-- Execute the monitors of the chain.
mcmcExecuteMonitors :: Algorithm a => a -> MCMC ()
mcmcExecuteMonitors a = do
  e <- ask
  let s = settings e
      vb = sVerbosity s
      t0 = startingTime e
      iTotal = burnInIterations (sBurnIn s) + fromNIterations (sNIterations s)
  mStdLog <- liftIO (aExecuteMonitors vb t0 iTotal a)
  forM_ mStdLog (mcmcOutB "   ")

mcmcIterate :: Algorithm a => Int -> a -> MCMC a
mcmcIterate n a
  | n < 0 = error "mcmcIterate: Number of iterations is negative."
  | n == 0 = return a
  | otherwise = do
    p <- sParallelizationMode . settings <$> ask
    a' <- liftIO $ aIterate p a
    mcmcExecuteMonitors a'
    mcmcIterate (n -1) a'

mcmcNewRun :: Algorithm a => a -> MCMC a
mcmcNewRun a = do
  s <- reader settings
  mcmcInfoB "Start new MCMC sampler."
  mcmcInfoB "Initial state."
  mcmcInfoB $ aStdMonitorHeader a
  mcmcExecuteMonitors a
  mcmcInfoB $ aSummarizeCycle a
  a' <- mcmcBurnIn a
  a'' <- mcmcResetAcceptance a'
  let i = fromNIterations $ sNIterations s
  mcmcInfoS $ "Run chain for " ++ show i ++ " iterations."
  mcmcInfoB $ aStdMonitorHeader a''
  mcmcIterate i a''

mcmcContinueRun :: Algorithm a => a -> MCMC a
mcmcContinueRun a = do
  s <- reader settings
  let iTotal = fromNIterations (sNIterations s) + burnInIterations (sBurnIn s)
  mcmcInfoB "Continuation of MCMC sampler."
  let iCurrent = aIteration a
  mcmcInfoS $ "Current iteration: " ++ show iCurrent ++ "."
  mcmcInfoS $ "Total iterations: " ++ show iTotal ++ "."
  let di = iTotal - iCurrent
  mcmcInfoB $ aSummarizeCycle a
  mcmcInfoS $ "Run chain for " ++ show di ++ " iterations."
  mcmcInfoB $ aStdMonitorHeader a
  mcmcIterate di a

mcmcBurnIn :: Algorithm a => a -> MCMC a
mcmcBurnIn a = do
  s <- reader settings
  case sBurnIn s of
    NoBurnIn -> do
      mcmcInfoS "No burn in."
      return a
    BurnInWithoutAutoTuning n -> do
      mcmcInfoS $ "Burn in for " <> show n <> " iterations."
      mcmcInfoS "Auto tuning is disabled."
      mcmcInfoB $ aStdMonitorHeader a
      a' <- mcmcIterate n a
      mcmcInfoB $ aSummarizeCycle a'
      mcmcInfoB "Burn in finished."
      return a'
    BurnInWithAutoTuning n t -> do
      mcmcInfoS $ "Burn in for " ++ show n ++ " iterations."
      mcmcInfoS $ "Auto tuning is enabled with a period of " ++ show t ++ "."
      mcmcInfoB $ aStdMonitorHeader a
      a' <- mcmcBurnInWithAutoTuning n t a
      mcmcInfoB "Burn in finished."
      return a'

-- Auto tune the proposals.
mcmcAutotune :: Algorithm a => a -> MCMC a
mcmcAutotune a = do
  mcmcDebugB "Auto tune."
  return $ aAutoTune a

mcmcBurnInWithAutoTuning :: Algorithm a => Int -> Int -> a -> MCMC a
mcmcBurnInWithAutoTuning b t a
  | b > t = do
    a' <- mcmcResetAcceptance a
    a'' <- mcmcIterate t a'
    mcmcDebugB $ aSummarizeCycle a''
    a''' <- mcmcAutotune a''
    mcmcDebugB $ aStdMonitorHeader a''
    mcmcBurnInWithAutoTuning (b - t) t a'''
  | otherwise = do
    a' <- mcmcResetAcceptance a
    a'' <- mcmcIterate b a'
    mcmcInfoB $ aSummarizeCycle a''
    mcmcInfoS $ "Acceptance rates calculated over the last " <> show b <> " iterations."
    return a''

mcmcInitialize :: Algorithm a => a -> MCMC a
mcmcInitialize a = do
  mcmcInfoS $ aName a ++ " algorithm."
  s <- settings <$> ask
  mcmcDebugB "Opening monitors."
  a' <- liftIO $ aOpenMonitors (sAnalysisName s) (sExecutionMode s) a
  mcmcDebugB "Monitors opened."
  return a'

-- Save the MCMC run.
mcmcSave :: Algorithm a => a -> MCMC ()
mcmcSave a = do
  s <- reader settings
  case sSaveMode s of
    NoSave -> mcmcInfoB "Do not save the Markov chain."
    SaveWithTrace n -> do
      mcmcInfoB "Save settings."
      liftIO $ settingsSave s
      let nm = sAnalysisName s
      mcmcInfoS $
        "Save compressed Markov chain with trace of length " ++ show n ++ "."
      mcmcInfoB "For long traces, or complex objects, this may take a while."
      liftIO $ aSave nm a
      mcmcInfoB "Markov chain saved."

-- Report and finish up.
mcmcClose :: Algorithm a => a -> MCMC a
mcmcClose a = do
  mcmcDebugB "Closing MCMC run."
  mcmcInfoB $ aSummarizeCycle a
  mcmcInfoS $ aName a ++ " algorithm finished."
  mcmcSave a
  ti <- reader startingTime
  te <- liftIO getCurrentTime
  let dt = te `diffUTCTime` ti
  mcmcInfoB $ "Wall clock run time: " <> renderDuration dt <> "."
  mcmcInfoS $ "End time: " <> renderTime te
  a' <- liftIO $ aCloseMonitors a
  h <- reader logHandle
  liftIO $ forM_ h hClose
  return a'

-- Initialize the run, execute the run, and close the run.
mcmcRun :: Algorithm a => a -> MCMC a
mcmcRun a = do
  mcmcDebugB "The settings are:"
  reader settings >>= mcmcDebugS . ppShow

  -- Initialize.
  a' <- mcmcInitialize a
  mcmcReportTime

  -- Execute.
  a'' <- mcmcExecute a'

  -- Close.
  mcmcClose a''

-- | Run an MCMC algorithm with given settings.
mcmc :: Algorithm a => Settings -> a -> IO a
mcmc s a = do
  settingsCheck s $ aIteration a
  e <- initializeEnvironment s
  runReaderT (mcmcRun a) e

-- | Continue an MCMC algorithm for the given number of iterations.
mcmcContinue :: Algorithm a => Int -> Settings -> a -> IO a
mcmcContinue dn s = mcmc s'
  where
    n' = NIterations $ fromNIterations (sNIterations s) + dn
    s' = s {sNIterations = n', sExecutionMode = Continue}
