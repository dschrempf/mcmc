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
    mcmcOutT,
    mcmcOutS,
    mcmcWarnT,
    mcmcWarnS,
    mcmcInfoT,
    mcmcInfoS,
    mcmcDebugT,
    mcmcDebugS,
    mcmcAutotune,
    mcmcResetA,
    mcmcSummarizeCycle,
    mcmcReport,
    mcmcMonitorStdOutHeader,
    mcmcMonitorExec,
    mcmcRun,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State
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

msgPrepare :: Char -> Text -> Text
msgPrepare c t = T.cons c $ ": " <> t

-- | Write to standard output and log file.
mcmcOutT :: Text -> Mcmc a ()
mcmcOutT msg = do
  h <- fromMaybe (error "mcmcOut: Log handle is missing.") <$> gets logHandle
  liftIO $ T.putStrLn msg >> T.hPutStrLn h msg

-- | Write to standard output and log file.
mcmcOutS :: String -> Mcmc a ()
mcmcOutS = mcmcOutT . T.pack

-- Perform warning action.
mcmcWarnA :: Mcmc a () -> Mcmc a ()
mcmcWarnA a = gets verbosity >>= \v -> info v a

-- | Print warning message.
mcmcWarnT :: Text -> Mcmc a ()
mcmcWarnT = mcmcWarnA . mcmcOutT . msgPrepare 'W'

-- | Print warning message.
mcmcWarnS :: String -> Mcmc a ()
mcmcWarnS = mcmcWarnT . T.pack

-- Perform info action.
mcmcInfoA :: Mcmc a () -> Mcmc a ()
mcmcInfoA a = gets verbosity >>= \v -> info v a

-- | Print info message.
mcmcInfoT :: Text -> Mcmc a ()
mcmcInfoT = mcmcInfoA . mcmcOutT . msgPrepare 'I'

-- | Print info message.
mcmcInfoS :: String -> Mcmc a ()
mcmcInfoS = mcmcInfoT . T.pack

-- Perform debug action.
mcmcDebugA :: Mcmc a () -> Mcmc a ()
mcmcDebugA a = gets verbosity >>= \v -> debug v a

-- | Print debug message.
mcmcDebugT :: Text -> Mcmc a ()
mcmcDebugT = mcmcDebugA . mcmcOutT . msgPrepare 'D'

-- | Print debug message.
mcmcDebugS :: String -> Mcmc a ()
mcmcDebugS = mcmcDebugT . T.pack

-- | Auto tune the 'Proposal's in the 'Cycle' of the chain. Reset acceptance counts.
-- See 'autotuneCycle'.
mcmcAutotune :: Mcmc a ()
mcmcAutotune = do
  mcmcDebugT "Auto tune."
  s <- get
  let a = acceptance s
      c = cycle s
      c' = autotuneCycle a c
  put $ s {cycle = c'}

-- | Reset acceptance counts.
mcmcResetA :: Mcmc a ()
mcmcResetA = do
  mcmcDebugT "Reset acceptance ratios."
  s <- get
  let a = acceptance s
  put $ s {acceptance = resetA a}

-- | Print short summary of 'Proposal's in 'Cycle'. See 'summarizeCycle'.
mcmcSummarizeCycle :: Mcmc a Text
mcmcSummarizeCycle = do
  a <- gets acceptance
  c <- gets cycle
  return $ summarizeCycle a c

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

-- Open log file.
mcmcOpenLog :: Mcmc a ()
mcmcOpenLog = do
  s <- get
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
  put s {logHandle = mh}
  mcmcDebugS $ "Log file name: " ++ lfn ++ "."
  mcmcDebugT "Log file opened."

-- | Set the total number of iterations, the current time and open the
-- 'Monitor's of the chain. See 'mOpen'.
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
      frc = forceOverwrite s
  m' <-
    if n == 0
      then liftIO $ mOpen nm frc m
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
  mcmcInfoT "Initial state."
  mcmcMonitorStdOutHeader
  mcmcMonitorExec

-- | Print header line of standard output monitor.
mcmcMonitorStdOutHeader :: Mcmc a ()
mcmcMonitorStdOutHeader = do
  m <- gets monitor
  mcmcInfoA $ mcmcOutT $ mHeader m

-- Save the status of an MCMC run. See 'saveStatus'.
mcmcSave :: ToJSON a => Mcmc a ()
mcmcSave = do
  s <- get
  if save s
    then do
      mcmcInfoT "Save Markov chain. For long chains, this may take a while."
      liftIO $ saveStatus (name s <> ".mcmc") s
      mcmcInfoT "Done saving Markov chain."
    else mcmcInfoT "Do not save the Markov chain."

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: ToJSON a => Mcmc a ()
mcmcMonitorExec = do
  s <- get
  let i = iteration s
      j = iterations s + fromMaybe 0 (burnInIterations s)
      m = monitor s
      (ss, st) = fromMaybe (error "mcmcMonitorExec: Starting state and time not set.") (start s)
      tr = trace s
      vb = verbosity s
  mt <- liftIO $ mExec vb i ss st tr j m
  forM_ mt mcmcOutT

-- | Close the 'Monitor's of the chain. See 'mClose'.
mcmcClose :: ToJSON a => Mcmc a ()
mcmcClose = do
  s <- get
  mcmcSummarizeCycle >>= mcmcInfoT
  mcmcInfoT "Metropolis-Hastings sampler finished."
  let m = monitor s
  m' <- liftIO $ mClose m
  put $ s {monitor = m'}
  mcmcSave
  t <- liftIO getCurrentTime
  let rt = case start s of
        Nothing -> error "mcmcClose: Start time not set."
        Just (_, st) -> t `diffUTCTime` st
  mcmcInfoT $ "Wall clock run time: " <> renderDuration rt <> "."
  mcmcInfoS $ "End time: " <> fTime t
  case logHandle s of
    Just h -> liftIO $ hClose h
    Nothing -> return ()

-- | Run an MCMC algorithm.
mcmcRun :: ToJSON a => Mcmc a () -> Status a -> IO (Status a)
mcmcRun algorithm = execStateT $ do mcmcInit
                                    algorithm
                                    mcmcClose
