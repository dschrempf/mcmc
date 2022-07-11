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
import Control.Monad.Trans.RWS.CPS
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
type MCMC a = RWST (Environment Settings) () a IO

mcmcExecute :: Algorithm a => MCMC a ()
mcmcExecute = do
  logDebugB "Executing MCMC run."
  s <- asks settings
  case sExecutionMode s of
    Fail -> mcmcNewRun
    Overwrite -> mcmcNewRun
    Continue -> mcmcContinueRun
  logDebugB "Executed MCMC run."

mcmcResetAcceptance :: Algorithm a => MCMC a ()
mcmcResetAcceptance = do
  logDebugB "Reset acceptance rates."
  modify aResetAcceptance

mcmcExceptionHandler :: Algorithm a => Environment Settings -> a -> AsyncException -> IO b
mcmcExceptionHandler e a err = do
  putStrLn ""
  putStrLn "INTERRUPT!"
  putStrLn "Try to terminate gracefully and save chain for continuation."
  putStrLn "Press CTRL-C (again) to terminate now."
  putStrLn "Closing output files."
  _ <- aCloseMonitors a
  closeEnvironment e
  putStrLn "Saving settings."
  let s = settings e
  settingsSave s
  putStrLn "Saving compressed MCMC analysis."
  putStrLn "For long traces, or complex objects, this may take a while."
  let nm = sAnalysisName s
  aSave nm a
  putStrLn "Markov chain saved. Analysis can be continued."
  putStrLn "Graceful termination successful."
  putStrLn "Rethrowing error."
  throw err

-- XXX: Exception handling. Is it enough to mask execution of monitors and catch
-- UserInterrupt during iterations?

mcmcExecuteMonitors :: Algorithm a => MCMC a ()
mcmcExecuteMonitors = do
  e <- ask
  let s = settings e
      vb = sVerbosity s
      t0 = startingTime e
      iTotal = burnInIterations (sBurnIn s) + fromIterations (sIterations s)
  a <- get
  -- NOTE: Mask asynchronous exceptions when writing monitor files.
  mStdLog <- liftIO $ mask_ $ aExecuteMonitors vb t0 iTotal a
  forM_ mStdLog (logOutB "   ")

mcmcIterate :: Algorithm a => IterationMode -> Int -> MCMC a ()
mcmcIterate m n
  | n < 0 = error "mcmcIterate: Number of iterations is negative."
  | n == 0 = pure ()
  | otherwise = do
      e <- ask
      p <- sParallelizationMode . settings <$> ask
      a <- get
      -- NOTE: User interrupt is handled during iterations.
      a' <- liftIO $ catch (aIterate m p a) (mcmcExceptionHandler e a)
      put a'
      mcmcExecuteMonitors
      mcmcIterate m (n - 1)

mcmcNewRun :: Algorithm a => MCMC a ()
mcmcNewRun = do
  s <- reader settings
  logInfoB "Start new MCMC sampler."
  logInfoB "Initial state."
  get >>= logInfoB . aStdMonitorHeader
  mcmcExecuteMonitors
  isInvalid <- gets aIsInvalidState
  when isInvalid (logWarnB "The initial state is invalid!")
  mcmcBurnIn
  logInfoS $ "Clean chain after burn in."
  let tl = sTraceLength s
  a <- get
  a' <- liftIO $ aCleanAfterBurnIn tl a
  put a'
  let i = fromIterations $ sIterations s
  logInfoS $ "Run chain for " ++ show i ++ " iterations."
  get >>= logInfoB . aStdMonitorHeader
  mcmcIterate AllProposals i

mcmcContinueRun :: Algorithm a => MCMC a ()
mcmcContinueRun = do
  s <- reader settings
  let iBurnIn = burnInIterations (sBurnIn s)
      iNormal = fromIterations (sIterations s)
      iTotal = iBurnIn + iNormal
  logInfoB "Continuation of MCMC sampler."
  iCurrent <- gets aIteration
  logInfoS $ "Burn in iterations: " ++ show iBurnIn ++ "."
  logInfoS $ "Normal iterations: " ++ show iNormal ++ "."
  logInfoS $ "Total iterations: " ++ show iTotal ++ "."
  logInfoS $ "Current iteration: " ++ show iCurrent ++ "."
  when (iCurrent < iBurnIn) $ error "mcmcContinueRun: Can not continue burn in."
  let di = iTotal - iCurrent
  get >>= logInfoB . aSummarizeCycle AllProposals
  logInfoS $ "Run chain for " ++ show di ++ " iterations."
  get >>= logInfoB . aStdMonitorHeader
  mcmcIterate AllProposals di

mcmcBurnIn :: Algorithm a => MCMC a ()
mcmcBurnIn = do
  s <- reader settings
  case sBurnIn s of
    NoBurnIn -> do
      get >>= logInfoB . aSummarizeCycle AllProposals
      logInfoS "No burn in."
    BurnInWithoutAutoTuning n -> do
      get >>= logInfoB . aSummarizeCycle AllProposals
      logInfoS $ "Burn in for " <> show n <> " iterations."
      logInfoS "Auto tuning is disabled."
      get >>= logInfoB . aStdMonitorHeader
      mcmcIterate AllProposals n
      get >>= logInfoB . aSummarizeCycle AllProposals
      mcmcResetAcceptance
      logInfoB "Burn in finished."
    BurnInWithAutoTuning n t -> do
      get >>= logInfoB . aSummarizeCycle AllProposals
      logInfoS $ "Burn in for " ++ show n ++ " iterations."
      logInfoS $ "Auto tuning is enabled with a period of " ++ show t ++ "."
      get >>= logInfoB . aStdMonitorHeader
      let (m, r) = n `divMod` t
          -- Don't add another auto tune period if r == 0, because then we auto
          -- tune without acceptance counts and get NaNs.
          xs = replicate m t <> [r | r > 0]
      mcmcBurnInWithAutoTuning AllProposals xs
      logInfoB "Burn in finished."
    BurnInWithCustomAutoTuning xs ys -> do
      logInfoS $ "Burn in for " ++ show (sum xs + sum ys) ++ " iterations."
      if null xs
        then get >>= logInfoB . aSummarizeCycle AllProposals
        else do
          get >>= logInfoB . aSummarizeCycle FastProposals
          logInfoS $ "Fast custom auto tuning with periods " ++ show xs ++ "."
          get >>= logInfoB . aStdMonitorHeader
          mcmcBurnInWithAutoTuning FastProposals xs
      logInfoS $ "Full custom auto tuning with periods " ++ show ys ++ "."
      get >>= logInfoB . aStdMonitorHeader
      mcmcBurnInWithAutoTuning AllProposals ys
      logInfoB "Burn in finished."

-- Auto tune the proposals.
mcmcAutotune :: Algorithm a => TuningType -> Int -> MCMC a ()
mcmcAutotune NormalTuningStep n = do
  logDebugB "Auto tune."
  a <- get
  a' <- liftIO $ aAutoTune NormalTuningStep n a
  put a'
mcmcAutotune LastTuningStep n = do
  logDebugB "Last auto tune."
  a <- get
  a' <- liftIO $ aAutoTune LastTuningStep n a
  put a'

mcmcBurnInWithAutoTuning :: Algorithm a => IterationMode -> [Int] -> MCMC a ()
mcmcBurnInWithAutoTuning _ [] = error "mcmcBurnInWithAutoTuning: Empty list."
mcmcBurnInWithAutoTuning m [x] = do
  -- Last round.
  mcmcIterate m x
  mcmcAutotune LastTuningStep x
  get >>= logInfoB . aSummarizeCycle m
  logInfoS $ "Acceptance rates calculated over the last " <> show x <> " iterations."
  mcmcResetAcceptance
mcmcBurnInWithAutoTuning m (x : xs) = do
  mcmcIterate m x
  mcmcAutotune NormalTuningStep x
  get >>= logDebugB . aSummarizeCycle m
  logDebugS $ "Acceptance rates calculated over the last " <> show x <> " iterations."
  get >>= logDebugB . aStdMonitorHeader
  mcmcResetAcceptance
  mcmcBurnInWithAutoTuning m xs

mcmcInitialize :: Algorithm a => MCMC a ()
mcmcInitialize = do
  a <- get
  logInfoS $ aName a ++ " algorithm."
  s <- settings <$> ask
  logDebugB "Opening monitors."
  a' <- liftIO $ aOpenMonitors (sAnalysisName s) (sExecutionMode s) a
  put a'
  logDebugB "Monitors opened."

-- Save the MCMC run.
mcmcSave :: Algorithm a => MCMC a ()
mcmcSave = do
  s <- reader settings
  case sSaveMode s of
    NoSave -> logInfoB "Do not save the MCMC analysis."
    Save -> do
      logInfoB "Save settings."
      liftIO $ settingsSave s
      let nm = sAnalysisName s
      logInfoB "Save compressed MCMC analysis."
      logInfoB "For long traces, or complex objects, this may take a while."
      a <- get
      liftIO $ aSave nm a
      logInfoB "Markov chain saved."

-- Report and finish up.
mcmcClose :: Algorithm a => MCMC a ()
mcmcClose = do
  logDebugB "Closing MCMC run."
  get >>= logInfoB . aSummarizeCycle AllProposals
  get >>= logInfoS . (\a -> aName a ++ " algorithm finished.")
  mcmcSave
  logInfoEndTime
  a <- get
  a' <- liftIO $ aCloseMonitors a
  put a'
  e <- ask
  liftIO $ closeEnvironment e

-- Initialize the run, execute the run, and close the run.
mcmcRun :: Algorithm a => MCMC a ()
mcmcRun = do
  -- Header.
  logInfoHeader
  reader settings >>= logInfoB . settingsPrettyPrint

  -- Initialize.
  mcmcInitialize
  logInfoStartingTime

  -- Execute.
  mcmcExecute

  -- Close.
  mcmcClose

-- | Run an MCMC algorithm with given settings.
mcmc :: Algorithm a => Settings -> a -> IO a
mcmc s a = do
  settingsCheck s $ aIteration a
  e <- initializeEnvironment s
  sn <$> runRWST mcmcRun e a
  where
    sn (_, x, _) = x

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
