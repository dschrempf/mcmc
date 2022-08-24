{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Mcmc
-- Description :  Framework for running Markov chain Monte Carlo samplers
-- Copyright   :  2021 Dominik Schrempf
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

import Control.Exception
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Reader
import Mcmc.Algorithm
import Mcmc.Cycle
import Mcmc.Environment
import Mcmc.Logger
import Mcmc.Proposal (TuningType (LastTuningStep, NormalTuningStep))
import Mcmc.Settings
import System.IO
import Prelude hiding (cycle)

-- The MCMC algorithm has read access to an environment and uses an algorithm
-- transforming the state @a@.
type MCMC = ReaderT (Environment Settings) IO

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

mcmcExceptionHandler :: Algorithm a => Environment Settings -> a -> AsyncException -> IO b
mcmcExceptionHandler e a err = do
  putStrLn ""
  putStrLn "INTERRUPT!"
  putStrLn "Try to terminate gracefully and save chain for continuation."
  putStrLn "Press CTRL-C (again) to terminate now."
  putStrLn "Closing output files."
  _ <- aCloseMonitors a
  closeEnvironment e
  runReaderT (mcmcSave a) e
  putStrLn "Graceful termination successful."
  putStrLn "Rethrowing error."
  throw err

mcmcExecuteMonitors :: Algorithm a => a -> MCMC ()
mcmcExecuteMonitors a = do
  e <- ask
  let s = settings e
      vb = sVerbosity s
      t0 = startingTime e
      iTotal = burnInIterations (sBurnIn s) + fromIterations (sIterations s)
  -- NOTE: Mask asynchronous exceptions when writing monitor files.
  mStdLog <- liftIO $ mask_ $ aExecuteMonitors vb t0 iTotal a
  forM_ mStdLog (logOutB "   ")

mcmcIterate :: Algorithm a => IterationMode -> Int -> a -> MCMC a
mcmcIterate m n a
  | n < 0 = error "mcmcIterate: Number of iterations is negative."
  | n == 0 = return a
  | otherwise = do
      e <- ask
      p <- sParallelizationMode . settings <$> ask
      -- NOTE: User interrupt is only handled during iterations.
      a' <- liftIO $ catch (aIterate m p a) (mcmcExceptionHandler e a)
      -- NOTE: We may want to mask execution of monitors and catch exceptions
      -- even then?
      mcmcExecuteMonitors a'
      mcmcIterate m (n - 1) a'

mcmcNewRun :: Algorithm a => a -> MCMC a
mcmcNewRun a = do
  s <- reader settings
  logInfoB "Starting new MCMC sampler."
  logInfoB "Initial state."
  logInfoB $ aStdMonitorHeader a
  mcmcExecuteMonitors a
  when (aIsInvalidState a) (logWarnB "The initial state is invalid!")
  a' <- mcmcBurnIn a
  logInfoS $ "Cleaning chain after burn in."
  let tl = sTraceLength s
  a'' <- liftIO $ aCleanAfterBurnIn tl a'
  logInfoS $ "Saving chain after burn in."
  mcmcSave a''
  let i = fromIterations $ sIterations s
  logInfoS $ "Running chain for " ++ show i ++ " iterations."
  logInfoB $ aStdMonitorHeader a''
  mcmcIterate AllProposals i a''

mcmcContinueRun :: Algorithm a => a -> MCMC a
mcmcContinueRun a = do
  s <- reader settings
  let iBurnIn = burnInIterations (sBurnIn s)
      iNormal = fromIterations (sIterations s)
      iTotal = iBurnIn + iNormal
  logInfoB "Continuation of MCMC sampler."
  let iCurrent = aIteration a
  logInfoS $ "Burn in iterations: " ++ show iBurnIn ++ "."
  logInfoS $ "Normal iterations: " ++ show iNormal ++ "."
  logInfoS $ "Total iterations: " ++ show iTotal ++ "."
  logInfoS $ "Current iteration: " ++ show iCurrent ++ "."
  when (iCurrent < iBurnIn) $ error "mcmcContinueRun: Can not continue burn in."
  let di = iTotal - iCurrent
  logInfoB $ aSummarizeCycle AllProposals a
  logInfoS $ "Running chain for " ++ show di ++ " iterations."
  logInfoB $ aStdMonitorHeader a
  mcmcIterate AllProposals di a

mcmcBurnIn :: Algorithm a => a -> MCMC a
mcmcBurnIn a = do
  s <- reader settings
  case sBurnIn s of
    NoBurnIn -> do
      logInfoB $ aSummarizeCycle AllProposals a
      logInfoS "No burn in."
      return a
    BurnInWithoutAutoTuning n -> do
      logInfoB $ aSummarizeCycle AllProposals a
      logInfoS $ "Burning in for " <> show n <> " iterations."
      logInfoS "Auto tuning is disabled."
      logInfoB $ aStdMonitorHeader a
      a' <- mcmcIterate AllProposals n a
      logInfoB $ aSummarizeCycle AllProposals a'
      a'' <- mcmcResetAcceptance a'
      logInfoB "Burn in finished."
      return a''
    BurnInWithAutoTuning n t -> do
      logInfoB $ aSummarizeCycle AllProposals a
      logInfoS $ "Burning in for " ++ show n ++ " iterations."
      logInfoS $ "Auto tuning is enabled with a period of " ++ show t ++ "."
      logInfoB $ aStdMonitorHeader a
      let (m, r) = n `divMod` t
          -- Don't add another auto tune period if r == 0, because then we auto
          -- tune without acceptance counts and get NaNs.
          xs = replicate m t <> [r | r > 0]
      a' <- mcmcBurnInWithAutoTuning AllProposals xs a
      logInfoB "Burn in finished."
      return a'
    BurnInWithCustomAutoTuning xs ys -> do
      logInfoS $ "Burning in for " ++ show (sum xs + sum ys) ++ " iterations."
      a' <-
        if null xs
          then do
            logInfoB $ aSummarizeCycle AllProposals a
            pure a
          else do
            logInfoB $ aSummarizeCycle FastProposals a
            logInfoS $ "Fast custom auto tuning with periods " ++ show xs ++ "."
            logInfoB $ aStdMonitorHeader a
            mcmcBurnInWithAutoTuning FastProposals xs a
      logInfoS $ "Full custom auto tuning with periods " ++ show ys ++ "."
      logInfoB $ aStdMonitorHeader a
      a'' <- mcmcBurnInWithAutoTuning AllProposals ys a'
      logInfoB "Burn in finished."
      return a''

-- Auto tune the proposals.
mcmcAutotune :: Algorithm a => TuningType -> Int -> a -> MCMC a
mcmcAutotune NormalTuningStep n a = do
  logDebugB "Intermediate auto tune."
  liftIO $ aAutoTune NormalTuningStep n a
mcmcAutotune LastTuningStep n a = do
  logDebugB "Last auto tune."
  liftIO $ aAutoTune LastTuningStep n a

mcmcBurnInWithAutoTuning :: Algorithm a => IterationMode -> [Int] -> a -> MCMC a
mcmcBurnInWithAutoTuning _ [] _ = error "mcmcBurnInWithAutoTuning: Empty list."
mcmcBurnInWithAutoTuning m [x] a = do
  -- Last round.
  a' <- mcmcIterate m x a
  a'' <- mcmcAutotune LastTuningStep x a'
  logInfoB $ aSummarizeCycle m a''
  logInfoS $ "Acceptance rates calculated over the last " <> show x <> " iterations."
  mcmcResetAcceptance a''
mcmcBurnInWithAutoTuning m (x : xs) a = do
  a' <- mcmcIterate m x a
  a'' <- mcmcAutotune NormalTuningStep x a'
  logDebugB $ aSummarizeCycle m a''
  logDebugS $ "Acceptance rates calculated over the last " <> show x <> " iterations."
  logDebugB $ aStdMonitorHeader a''
  a''' <- mcmcResetAcceptance a''
  mcmcBurnInWithAutoTuning m xs a'''

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
    NoSave -> logInfoB "NoSave set; Do not save the MCMC analysis."
    Save -> do
      logInfoB "Saving settings."
      liftIO $ settingsSave s
      let nm = sAnalysisName s
      logInfoB "Saving compressed MCMC analysis."
      logInfoB "For long traces, or complex objects, this may take a while."
      liftIO $ aSave nm a
      logInfoB "Markov chain saved. Analysis can be continued."

-- Report and finish up.
mcmcClose :: Algorithm a => a -> MCMC a
mcmcClose a = do
  logDebugB "Closing MCMC run."
  logInfoB $ aSummarizeCycle AllProposals a
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
  -- Header.
  logInfoHeader
  reader settings >>= logInfoB . settingsPrettyPrint

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
-- process is killed or fails.
--
-- See:
--
-- - 'Mcmc.Algorithm.MHG.mhgLoad'
--
-- - 'Mcmc.Algorithm.MC3.mc3Load'
mcmcContinue :: Algorithm a => Iterations -> Settings -> a -> IO a
mcmcContinue dn s = mcmc s'
  where
    n' = Iterations $ fromIterations (sIterations s) + fromIterations dn
    s' = s {sIterations = n', sExecutionMode = Continue}
