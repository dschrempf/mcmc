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
    mcmcOut,
    mcmcWarn,
    mcmcWarnS,
    mcmcInfoClean,
    mcmcInfo,
    mcmcInfoS,
    mcmcDebug,
    mcmcDebugS,
    mcmcAutotune,
    mcmcResetA,
    mcmcSummarizeCycle,
    mcmcInit,
    mcmcReport,
    mcmcMonitorStdOutHeader,
    mcmcMonitorExec,
    mcmcClose,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict
import Data.Aeson
import Data.Maybe
import qualified Data.Text.Lazy as T
import Data.Text.Lazy (Text)
import qualified Data.Text.Lazy.IO as T
import Data.Time.Clock
import Data.Time.Format
import Mcmc.Monitor
import Mcmc.Monitor.Time
import Mcmc.Proposal
import Mcmc.Save
import Mcmc.Status hiding (debug)
import Mcmc.Verbosity
import System.Directory
import System.IO
import Prelude hiding (cycle)

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO

-- | Write to standard output and log file.
mcmcOut :: Text -> Mcmc a ()
mcmcOut msg = do
  h <- fromMaybe (error "mcmcOut: Log handle is missing.") <$> gets logHandle
  liftIO $ T.putStrLn msg >> T.hPutStrLn h msg

-- | Print warning.
mcmcWarn :: Text -> Mcmc a ()
mcmcWarn msg = do
  v <- gets verbosity
  let msg' = T.intercalate "\n" $ map ("-- WARNING:" <>) $ T.lines msg
  warn v $ mcmcOut msg'

-- | Print warning.
mcmcWarnS :: String -> Mcmc a ()
mcmcWarnS = mcmcWarn . T.pack

-- | Print info message without line prefix.
mcmcInfoClean :: Text -> Mcmc a ()
mcmcInfoClean msg = do
  v <- gets verbosity
  info v $ mcmcOut msg

-- | Print info message.
mcmcInfo :: Text -> Mcmc a ()
mcmcInfo = mcmcInfoClean . T.intercalate "\n" . map ("-- " <>) . T.lines

-- | Print info message.
mcmcInfoS :: String -> Mcmc a ()
mcmcInfoS = mcmcInfo . T.pack

-- | Print debug message.
mcmcDebug :: Text -> Mcmc a ()
mcmcDebug msg = do
  v <- gets verbosity
  let msg' = T.intercalate "\n" $ map ("-- DEBUG: " <>) $ T.lines msg
  debug v $ mcmcOut msg'

-- | Print debug message.
mcmcDebugS :: String -> Mcmc a ()
mcmcDebugS = mcmcDebug . T.pack

-- | Auto tune the 'Proposal's in the 'Cycle' of the chain. Reset acceptance counts.
-- See 'autotuneCycle'.
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

-- | Print short summary of 'Proposal's in 'Cycle'. If 'True', also print acceptance
-- ratios. See 'summarizeCycle'.
mcmcSummarizeCycle :: Mcmc a Text
mcmcSummarizeCycle = do
  a <- gets acceptance
  c <- gets cycle
  return $ summarizeCycle a c

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

-- | Set the total number of iterations, the current time and open the
-- 'Monitor's of the chain. See 'mOpen'.
mcmcInit :: Mcmc a ()
mcmcInit = do
  s <- get
  -- Log file.
  let lfn = name s ++ ".log"
      n = iteration s
      frc = forceOverwrite s
  fe <- liftIO $ doesFileExist lfn
  mh <- liftIO $ case verbosity s of
    Quiet -> return Nothing
    _ -> Just <$> case (fe, n, frc) of
      (False, _, _) -> openFile lfn WriteMode
      (True, 0, True) -> openFile lfn WriteMode
      (True, 0, False) -> error "mcmcInit: Log file exists; use 'force' to overwrite output files."
      (True, _, _) -> openFile lfn AppendMode
  let s' = s {logHandle = mh}
  put s'
  -- Start time.
  t <- liftIO getCurrentTime
  mcmcInfoS $ "Start time: " <> fTime t
  -- Monitor.
  let m = monitor s'
      nm = name s'
  m' <-
    if n == 0
      then liftIO $ mOpen nm frc m
      else liftIO $ mAppend nm m
  put $ s' {monitor = m', start = Just (n, t)}

-- | Report what is going to be done.
mcmcReport :: ToJSON a => Mcmc a ()
mcmcReport = do
  s <- get
  let b = burnInIterations s
      t = autoTuningPeriod s
      n = iterations s
  case b of
    Just b' -> mcmcInfoS $ "Burn in for " <> show b' <> " iterations."
    Nothing -> return ()
  case t of
    Just t' ->
      mcmcInfoS $
        "Auto tune every "
          <> show t'
          <> " iterations (during burn in only)."
    Nothing -> return ()
  mcmcInfoS $ "Run chain for " <> show n <> " iterations."
  mcmcInfo "Initial state."
  mcmcMonitorStdOutHeader
  mcmcMonitorExec

-- | Print header line of standard output monitor.
mcmcMonitorStdOutHeader :: Mcmc a ()
mcmcMonitorStdOutHeader = do
  m <- gets monitor
  v <- gets verbosity
  info v $ mcmcOut $ mHeader m

-- Save the status of an MCMC run. See 'saveStatus'.
mcmcSave :: ToJSON a => Mcmc a ()
mcmcSave = do
  s <- get
  if save s
    then do
      mcmcInfo "Save Markov chain. For long chains, this may take a while."
      liftIO $ saveStatus (name s <> ".mcmc") s
      mcmcInfo "Done saving Markov chain."
    else mcmcInfo "Do not save the Markov chain."

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
  _ <- mcmcInfoClean <$> mcmcSummarizeCycle
  mcmcInfo "Metropolis-Hastings sampler finished."
  let m = monitor s
  m' <- liftIO $ mClose m
  put $ s {monitor = m'}
  mcmcSave
  t <- liftIO getCurrentTime
  let rt = case start s of
        Nothing -> error "mcmcClose: Start time not set."
        Just (_, st) -> t `diffUTCTime` st
  mcmcInfo $ "Wall clock run time: " <> renderDuration rt <> "."
  mcmcInfoS $ "End time: " <> fTime t
  case logHandle s of
    Just h -> liftIO $ hClose h
    Nothing -> return ()
