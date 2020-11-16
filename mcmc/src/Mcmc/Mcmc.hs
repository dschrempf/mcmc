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
import Mcmc.Chain
import Mcmc.Environment
import Mcmc.Monitor
import Mcmc.Monitor.Time
import Mcmc.Proposal
import Mcmc.Save
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
  h <- fromMaybe (error "mcmcOut: Log handle is missing.") <$> reader envLogHandle
  liftIO $ BL.putStrLn msg >> BL.hPutStrLn h msg

-- | Write to standard output and log file.
mcmcOutS :: String -> Mcmc a ()
mcmcOutS = mcmcOutB . BL.pack

-- Perform warning action.
mcmcWarnA :: Mcmc a () -> Mcmc a ()
mcmcWarnA a = reader envVerbosity >>= \v -> when (v >= Warn) a

-- | Print warning message.
mcmcWarnB :: BL.ByteString -> Mcmc a ()
mcmcWarnB = mcmcWarnA . mcmcOutB . msgPrepare 'W'

-- | Print warning message.
mcmcWarnS :: String -> Mcmc a ()
mcmcWarnS = mcmcWarnB . BL.pack

-- Perform info action.
mcmcInfoA :: Mcmc a () -> Mcmc a ()
mcmcInfoA a = reader envVerbosity >>= \v -> when (v >= Info) a

-- | Print info message.
mcmcInfoB :: BL.ByteString -> Mcmc a ()
mcmcInfoB = mcmcInfoA . mcmcOutB . msgPrepare 'I'

-- | Print info message.
mcmcInfoS :: String -> Mcmc a ()
mcmcInfoS = mcmcInfoB . BL.pack

-- Perform debug action.
mcmcDebugA :: Mcmc a () -> Mcmc a ()
mcmcDebugA a = reader envVerbosity >>= \v -> when (v == Debug) a

-- | Print debug message.
mcmcDebugB :: BL.ByteString -> Mcmc a ()
mcmcDebugB = mcmcDebugA . mcmcOutB . msgPrepare 'D'

-- | Print debug message.
mcmcDebugS :: String -> Mcmc a ()
mcmcDebugS = mcmcDebugB . BL.pack

-- | Auto tune the 'Proposal's in the 'Cycle' of the chain. Reset acceptance counts.
-- See 'autotuneCycle'.
mcmcAutotune :: Algorithm a => Mcmc a ()
mcmcAutotune = do
  mcmcDebugB "Auto tune."
  modify autotune

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
mcmcResetA :: Mcmc a ()
mcmcResetA = do
  mcmcDebugB "Reset acceptance ratios."
  s <- get
  let a = acceptance s
  put $ s {acceptance = resetA a}

-- | Print short summary of 'Proposal's in 'Cycle'. See 'summarizeCycle'.
mcmcSummarizeCycle :: Mcmc a BL.ByteString
mcmcSummarizeCycle = do
  a <- gets acceptance
  c <- gets cycle
  return $ summarizeCycle a c

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

-- Set the total number of iterations, the current time and open the 'Monitor's
-- of the chain. See 'mOpen'.
mcmcInit :: Mcmc a ()
mcmcInit = do
  mcmcOpenLog
  s <- get
  -- Start time.
  t <- liftIO getCurrentTime
  mcmcInfoS $ "Start time: " <> fTime t
  -- Monitor.
  let m = monitor s
      n = iteration s
      nm = name s
  frc <- reader overwrite
  m' <- if n == 0 then liftIO $ mOpen nm frc m else liftIO $ mAppend nm m
  put $ s {monitor = m', start = Just (n, t)}

-- | Report what is going to be done.
mcmcReport :: ToJSON a => Mcmc a ()
mcmcReport = do
  s <- get
  let b = burnInIterations s
      t = autoTuningPeriod s
      n = iterations s
      c = cleaner s
  case b of
    Just b' -> mcmcInfoS $ "Burn in for " <> show b' <> " iterations."
    Nothing -> return ()
  case t of
    Just t' -> mcmcInfoS $ "Auto tune every " <> show t' <> " iterations (during burn in only)."
    Nothing -> return ()
  case c of
    Just (Cleaner c' _) -> mcmcInfoS $ "Clean state every " <> show c' <> " iterations."
    Nothing -> return ()
  mcmcInfoS $ "Run chain for " <> show n <> " iterations."
  mcmcInfoB "Initial state."
  mcmcMonitorExec

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

-- | Run an MCMC algorithm.
mcmcRun :: ToJSON a => Mcmc a () -> Environment -> Chain a -> IO (Chain a, ())
mcmcRun algorithm = execRWST $ do
  mcmcInit
  algorithm
  mcmcClose
