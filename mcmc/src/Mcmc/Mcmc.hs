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
-- Functions to work with the 'Mcmc' state transformer.
module Mcmc.Mcmc
  ( Mcmc,
    mcmcOutB,
    mcmcOutS,
    mcmcWarnB,
    mcmcWarnS,
    mcmcInfoB,
    mcmcInfoS,
    mcmcDebugB,
    mcmcDebugS,
    mcmcAutotune,
    mcmcClean,
    mcmcResetA,
    mcmcSummarizeCycle,
    mcmcReport,
    mcmcMonitorExec,
    mcmcRun,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.RWS.CPS
import Data.Aeson
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Maybe
import Data.Time.Clock
import Data.Time.Format
import Mcmc.Algorithm
import Mcmc.Environment
import Mcmc.Monitor
import Mcmc.Monitor.Time
import Mcmc.Save
import Mcmc.Settings
import Numeric.Log
import System.IO
import Prelude hiding (cycle)

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = RWST Environment () a IO

msgPrepare :: Char -> BL.ByteString -> BL.ByteString
msgPrepare c t = BL.cons c $ ": " <> t

-- | Write to standard output and log file.
mcmcOutB :: BL.ByteString -> Mcmc a ()
mcmcOutB msg = do
  h <- fromMaybe (error "mcmcOut: Log handle is missing.") <$> reader logHandle
  liftIO $ BL.putStrLn msg >> BL.hPutStrLn h msg

-- | Write to standard output and log file.
mcmcOutS :: String -> Mcmc a ()
mcmcOutS = mcmcOutB . BL.pack

-- Perform warning action.
mcmcWarnA :: Mcmc a () -> Mcmc a ()
mcmcWarnA a = reader verbosity >>= \v -> when (v >= Warn) a

-- | Print warning message.
mcmcWarnB :: BL.ByteString -> Mcmc a ()
mcmcWarnB = mcmcWarnA . mcmcOutB . msgPrepare 'W'

-- | Print warning message.
mcmcWarnS :: String -> Mcmc a ()
mcmcWarnS = mcmcWarnB . BL.pack

-- Perform info action.
mcmcInfoA :: Mcmc a () -> Mcmc a ()
mcmcInfoA a = reader verbosity >>= \v -> when (v >= Info) a

-- | Print info message.
mcmcInfoB :: BL.ByteString -> Mcmc a ()
mcmcInfoB = mcmcInfoA . mcmcOutB . msgPrepare 'I'

-- | Print info message.
mcmcInfoS :: String -> Mcmc a ()
mcmcInfoS = mcmcInfoB . BL.pack

-- Perform debug action.
mcmcDebugA :: Mcmc a () -> Mcmc a ()
mcmcDebugA a = reader verbosity >>= \v -> when (v == Debug) a

-- | Print debug message.
mcmcDebugB :: BL.ByteString -> Mcmc a ()
mcmcDebugB = mcmcDebugA . mcmcOutB . msgPrepare 'D'

-- | Print debug message.
mcmcDebugS :: String -> Mcmc a ()
mcmcDebugS = mcmcDebugB . BL.pack

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

-- Set the total number of iterations, the current time and open the 'Monitor's
-- of the chain. See 'mOpen'.
mcmcInit :: Algorithm a => Mcmc a ()
mcmcInit = do
  t <- reader startingTime
  mcmcInfoS $ "Starting time: " <> fTime t

-- Report what is going to be done.
mcmcStartingReport :: ToJSON a => Mcmc a ()
mcmcStartingReport = do
  s <- settings <$> ask
  let b = burnIn s
      i = iterations s
      c = cleaner s
  case b of
    NoBurnIn ->
      mcmcInfoS "No burn in."
    BurnInNoAutoTuning n ->
      mcmcInfoS $ "Burn in for " <> show n <> " iterations without auto tuning."
    BurnInWithAutoTuning n t ->
      mcmcInfoS $
        "Burn in for "
          <> show n
          <> " iterations with auto tuning every "
          <> show t " iterations."
  mcmcInfoS $ "Run chain for " <> show i <> " iterations."
  case c of
    NoClean -> return ()
    CleanEvery n -> mcmcInfoS $ "Clean state every " <> show n <> " iterations."
  mcmcInfoB "Initial state."
  mcmcMonitorExec

-- | Auto tune the 'Proposal's in the 'Cycle' of the chain. Reset acceptance counts.
-- See 'autoTuneCycle'.
mcmcAutotune :: Algorithm a => Mcmc a ()
mcmcAutotune = do
  mcmcDebugB "Auto tune."
  modify autoTune

-- | Clean the state.
mcmcClean :: Mcmc a ()
mcmcClean = do
  mcmcDebugB "Clean state."
  (pr, lh) <- report <$> get
  mcmcDebugS $
    "Current log prior and log likelihood: " ++ show (ln pr) ++ ", " ++ show (ln lh) ++ "."
  modify clean
  (pr', lh') <- report <$> get
  mcmcDebugS $
    "New log prior and log likelihood: " ++ show (ln pr') ++ ", " ++ show (ln lh') ++ "."
  let dLogPr = abs $ ln pr - ln pr'
      dLogLh = abs $ ln lh - ln lh'
  when
    (dLogPr > 1e-4)
    (mcmcWarnS $ "The logarithms of old and new log prior differ by " ++ show dLogPr ++ ".")
  when
    (dLogPr > 1e-4)
    (mcmcWarnS $ "The logarithms of old and new likelihood differ by " ++ show dLogLh ++ ".")

-- | Reset acceptance counts.
mcmcResetA :: Algorithm a => Mcmc a ()
mcmcResetA = do
  mcmcDebugB "Reset acceptance ratios."
  modify resetAcceptance

-- | Print short summary of 'Proposal's in 'Cycle'. See 'summarizeCycle'.
mcmcSummarizeCycle :: Algorithm a => Mcmc a BL.ByteString
mcmcSummarizeCycle = gets summarizeCycle

-- Save the status of an MCMC run. See 'saveStatus'.
mcmcSave :: ToJSON a => Mcmc a ()
mcmcSave = do
  env <- ask
  st <- get
  case saveMode env of
    NoSave -> mcmcInfoB "Do not save the Markov chain."
    SaveWithTrace n -> do
      mcmcInfoB $ "Save Markov chain with trace of length " <> BL.pack (show n) <> "."
      mcmcInfoB "For long traces, or complex objects, this may take a while."
      liftIO $ save env st
      mcmcInfoB "Done saving Markov chain."

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: ToJSON a => Mcmc a ()
mcmcMonitorExec = do
  vb <- reader verbosity
  s <- get
  let i = iteration s
      j = iterations s + fromMaybe 0 (burnInIterations s)
      m = monitor s
      (ss, st) = fromMaybe (error "mcmcMonitorExec: Starting state and time not set.") (start s)
      tr = trace s
  mt <- liftIO $ mExec vb i ss st tr j m
  forM_ mt mcmcOutB

-- Close the 'Monitor's of the chain. See 'mClose'.
mcmcClose :: ToJSON a => Mcmc a ()
mcmcClose = do
  s <- get
  mcmcSummarizeCycle >>= mcmcInfoB
  mcmcInfoB "Metropolis-Hastings sampler finished."
  let m = monitor s
  m' <- liftIO $ mClose m
  put $ s {monitor = m'}
  mcmcSave
  t <- liftIO getCurrentTime
  let rt = case start s of
        Nothing -> error "mcmcClose: Start time not set."
        Just (_, st) -> t `diffUTCTime` st
  mcmcInfoB $ "Wall clock run time: " <> renderDuration rt <> "."
  mcmcInfoS $ "End time: " <> fTime t
  case logHandle s of
    Just h -> liftIO $ hClose h
    Nothing -> return ()

mcmcRun :: Algorithm a => Mcmc a ()
mcmcRun = do
  mcmcInit
  mcmcStartingReport
  -- TODO!
  undefined
  mcmcClose

-- | Run an MCMC algorithm.
mcmcExec :: Algorithm a => Environment -> a -> IO (a, ())
mcmcExec env' alg = do
  env <- openLogFile
  -- TODO.
  -- mcmcDebugS $ "Log file name: " ++ lfn ++ "."
  -- mcmcDebugB "Log file opened."
  execRWST mcmcRun env alg
