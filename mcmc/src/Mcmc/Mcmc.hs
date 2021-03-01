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
import Mcmc.Algorithm
import Mcmc.Environment
import Mcmc.Logger
import Mcmc.Settings
import System.IO
import Text.Show.Pretty
import Prelude hiding (cycle)

-- The MCMC algorithm has read access to an environment and uses an algorithm
-- transforming the state @a@.
type MCMC a = ReaderT (Environment Settings) IO a

mcmcExecute :: Algorithm a => a -> MCMC a
mcmcExecute a = do
  logDebugB "Executing MCMC run."
  s <- reader settings
  a' <- case sExecutionMode s of
    Fail -> mcmcNewRun a
    Overwrite -> mcmcNewRun a
    Continue -> mcmcContinueRun a
  logDebugB "Executed MCMC run."
  return a'

mcmcResetAcceptance :: Algorithm a => a -> MCMC a
mcmcResetAcceptance a = do
  logDebugB "Reset acceptance rates."
  return $ aResetAcceptance a

mcmcExecuteMonitors :: Algorithm a => a -> MCMC ()
mcmcExecuteMonitors a = do
  e <- ask
  let s = settings e
      vb = sVerbosity s
      t0 = startingTime e
      iTotal = burnInIterations (sBurnIn s) + fromIterations (sIterations s)
  mStdLog <- liftIO (aExecuteMonitors vb t0 iTotal a)
  forM_ mStdLog (logOutB "   ")

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
  logInfoB "Start new MCMC sampler."
  logInfoB "Initial state."
  logInfoB $ aStdMonitorHeader a
  mcmcExecuteMonitors a
  logInfoB $ aSummarizeCycle a
  a' <- mcmcBurnIn a
  let i = fromIterations $ sIterations s
  logInfoS $ "Run chain for " ++ show i ++ " iterations."
  logInfoB $ aStdMonitorHeader a'
  mcmcIterate i a'

mcmcContinueRun :: Algorithm a => a -> MCMC a
mcmcContinueRun a = do
  s <- reader settings
  let iTotal = fromIterations (sIterations s) + burnInIterations (sBurnIn s)
  logInfoB "Continuation of MCMC sampler."
  let iCurrent = aIteration a
  logInfoS $ "Current iteration: " ++ show iCurrent ++ "."
  logInfoS $ "Total iterations: " ++ show iTotal ++ "."
  let di = iTotal - iCurrent
  logInfoB $ aSummarizeCycle a
  logInfoS $ "Run chain for " ++ show di ++ " iterations."
  logInfoB $ aStdMonitorHeader a
  mcmcIterate di a

mcmcBurnIn :: Algorithm a => a -> MCMC a
mcmcBurnIn a = do
  s <- reader settings
  case sBurnIn s of
    NoBurnIn -> do
      logInfoS "No burn in."
      return a
    BurnInWithoutAutoTuning n -> do
      logInfoS $ "Burn in for " <> show n <> " iterations."
      logInfoS "Auto tuning is disabled."
      logInfoB $ aStdMonitorHeader a
      a' <- mcmcIterate n a
      logInfoB $ aSummarizeCycle a'
      a'' <- mcmcResetAcceptance a'
      logInfoB "Burn in finished."
      return a''
    BurnInWithAutoTuning n t -> do
      logInfoS $ "Burn in for " ++ show n ++ " iterations."
      logInfoS $ "Auto tuning is enabled with a period of " ++ show t ++ "."
      logInfoB $ aStdMonitorHeader a
      let (m, r) = n `divMod` t
          xs = replicate m t <> [r]
      a' <- mcmcBurnInWithAutoTuning xs a
      logInfoB "Burn in finished."
      return a'
    BurnInWithCustomAutoTuning xs -> do
      logInfoS $ "Burn in for " ++ show (sum xs) ++ " iterations."
      logInfoS $ "Custom auto tuning is enabled with a periods " ++ show xs ++ "."
      logInfoB $ aStdMonitorHeader a
      a' <- mcmcBurnInWithAutoTuning xs a
      logInfoB "Burn in finished."
      return a'

-- Auto tune the proposals.
mcmcAutotune :: Algorithm a => a -> MCMC a
mcmcAutotune a = do
  logDebugB "Auto tune."
  return $ aAutoTune a

mcmcBurnInWithAutoTuning :: Algorithm a => [Int] -> a -> MCMC a
mcmcBurnInWithAutoTuning [] _ = error "mcmcBurnInWithAutoTuning: Empty lisst."
mcmcBurnInWithAutoTuning [x] a = do
  -- Last round.
  a' <- mcmcIterate x a
  a'' <- mcmcAutotune a'
  logInfoB $ aSummarizeCycle a''
  logInfoS $ "Acceptance rates calculated over the last " <> show x <> " iterations."
  mcmcResetAcceptance a
mcmcBurnInWithAutoTuning (x:xs) a = do
  a' <- mcmcIterate x a
  a'' <- mcmcAutotune a'
  logDebugB $ aSummarizeCycle a''
  logDebugS $ "Acceptance rates calculated over the last " <> show x <> " iterations."
  logDebugB $ aStdMonitorHeader a''
  a''' <- mcmcResetAcceptance a
  mcmcBurnInWithAutoTuning xs a'''

mcmcInitialize :: Algorithm a => a -> MCMC a
mcmcInitialize a = do
  logInfoS $ aName a ++ " algorithm."
  s <- settings <$> ask
  logDebugB "Opening monitors."
  a' <- liftIO $ aOpenMonitors (sAnalysisName s) (sExecutionMode s) a
  logDebugB "Monitors opened."
  return a'

-- Save the MCMC run.
mcmcSave :: Algorithm a => a -> MCMC ()
mcmcSave a = do
  s <- reader settings
  case sSaveMode s of
    NoSave -> logInfoB "Do not save the MCMC analysis."
    Save -> do
      logInfoB "Save settings."
      liftIO $ settingsSave s
      let nm = sAnalysisName s
      logInfoB "Save compressed MCMC analysis."
      logInfoB "For long traces, or complex objects, this may take a while."
      liftIO $ aSave nm a
      logInfoB "Markov chain saved."

-- Report and finish up.
mcmcClose :: Algorithm a => a -> MCMC a
mcmcClose a = do
  logDebugB "Closing MCMC run."
  logInfoB $ aSummarizeCycle a
  logInfoS $ aName a ++ " algorithm finished."
  mcmcSave a
  logInfoEndTime
  a' <- liftIO $ aCloseMonitors a
  e <- ask
  liftIO $ closeEnvironment e
  return a'

-- Initialize the run, execute the run, and close the run.
mcmcRun :: Algorithm a => a -> MCMC a
mcmcRun a = do
  logDebugB "The MCMC settings are:"
  reader settings >>= logDebugS . ppShow

  -- Initialize.
  a' <- mcmcInitialize a
  logInfoStartingTime

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
--
-- Currently, it is only possible to continue MCMC algorithms that have
-- completed successfully. This restriction is necessary, because for parallel
-- chains, it is hardly possible to ensure all chains are synchronized when the
-- process is killed.
--
-- See:
--
-- - 'Mcmc.Algorithm.MHG.mhgLoad'
--
-- - 'Mcmc.Algorithm.MC3.mc3Load'
mcmcContinue :: Algorithm a => Int -> Settings -> a -> IO a
mcmcContinue dn s = mcmc s'
  where
    n' = Iterations $ fromIterations (sIterations s) + dn
    s' = s {sIterations = n', sExecutionMode = Continue}
