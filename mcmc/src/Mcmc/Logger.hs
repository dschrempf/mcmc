{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Logger
-- Description :  Minimal monad logger
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Jan 12 09:03:04 2021.
module Mcmc.Logger
  ( LogMode (..),
    Verbosity (..),
    HasLock (..),
    HasLogHandles (..),
    HasStartingTime (..),
    HasLogMode (..),
    HasVerbosity (..),
    Logger,
    logOutB,
    logDebugB,
    logDebugS,
    logWarnB,
    logWarnS,
    logInfoB,
    logInfoS,
    logInfoHeader,
    logInfoStartingTime,
    logInfoEndTime,
  )
where

import Control.Concurrent.MVar
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Reader
import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Time.Clock
import Mcmc.Monitor.Time
import System.IO
import Paths_mcmc (version)
import Data.Version (showVersion)

-- | Define where the log output should be directed to.
--
-- Logging is disabled if 'Verbosity' is set to 'Quiet'.
data LogMode = LogStdOutAndFile | LogStdOutOnly | LogFileOnly
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''LogMode)

-- | Not much to say here.
data Verbosity = Quiet | Warn | Info | Debug
  deriving (Eq, Ord, Read, Show)

$(deriveJSON defaultOptions ''Verbosity)

-- | Types with an output lock for concurrent output.
class HasLock e where
  getLock :: e -> MVar ()

-- | Types with logging information.
class HasLogHandles e where
  getLogHandles :: e -> [Handle]

-- | Types with starting time.
class HasStartingTime s where
  getStartingTime :: s -> UTCTime

-- | Types with a log mode.
class HasLogMode s where
  getLogMode :: s -> LogMode

-- | Types with verbosity.
class HasVerbosity s where
  getVerbosity :: s -> Verbosity

-- | Reader transformer used for logging to a file and to standard output.
type Logger e a = ReaderT e IO a

msgPrepare :: BL.ByteString -> BL.ByteString -> BL.ByteString
msgPrepare pref msg = BL.intercalate "\n" $ map (BL.append pref) $ BL.lines msg

-- Make sure that concurrent output is not scrambled.
atomicAction :: HasLock e => IO () -> Logger e ()
atomicAction a = do
  l <- reader getLock
  liftIO $ withMVar l (const a)

-- | Write to standard output and maybe to log file.
logOutB ::
  (HasLogHandles e, HasLock e) =>
  -- | Prefix.
  BL.ByteString ->
  -- | Message.
  BL.ByteString ->
  Logger e ()
logOutB pref msg = do
  hs <- reader getLogHandles
  mapM_ (atomicAction . (`BL.hPutStrLn` msg')) hs
  where
    msg' = msgPrepare pref msg

-- Perform debug action.
logDebugA :: (HasLock e, HasLogHandles e, HasVerbosity e) => Logger e () -> Logger e ()
logDebugA a = reader getVerbosity >>= \v -> when (v >= Debug) a

-- | Log debug message.
logDebugB :: (HasLock e, HasLogHandles e, HasVerbosity e) => BL.ByteString -> Logger e ()
logDebugB = logDebugA . logOutB "D: "

-- | Log debug message.
logDebugS :: (HasLock e, HasLogHandles e, HasVerbosity e) => String -> Logger e ()
logDebugS = logDebugB . BL.pack

-- Perform warning action.
logWarnA :: (HasLogHandles e, HasVerbosity e) => Logger e () -> Logger e ()
logWarnA a = reader getVerbosity >>= \v -> when (v >= Warn) a

-- | Log warning message.
logWarnB :: (HasLock e, HasLogHandles e, HasVerbosity e) => BL.ByteString -> Logger e ()
logWarnB = logWarnA . logOutB "W: "

-- | Log warning message.
logWarnS :: (HasLock e, HasLogHandles e, HasVerbosity e) => String -> Logger e ()
logWarnS = logWarnB . BL.pack

-- Perform info action.
logInfoA :: (HasLogHandles e, HasVerbosity e) => Logger e () -> Logger e ()
logInfoA a = reader getVerbosity >>= \v -> when (v >= Info) a

-- | Log info message.
logInfoB :: (HasLock e, HasLogHandles e, HasVerbosity e) => BL.ByteString -> Logger e ()
logInfoB = logInfoA . logOutB "   "

-- | Log info message.
logInfoS :: (HasLock e, HasLogHandles e, HasVerbosity e) => String -> Logger e ()
logInfoS = logInfoB . BL.pack

-- | Log info header.
logInfoHeader :: (HasLock e, HasLogHandles e, HasVerbosity e) => Logger e ()
logInfoHeader = do
  logInfoS (replicate 70 '-')
  logInfoS ("MCMC sampler; version " ++ showVersion version <> ".")
  logInfoS "Developed by: Dominik Schrempf."
  logInfoS "License: GPL-3.0-or-later."
  logInfoS (replicate 70 '-')

-- | Log starting time.
logInfoStartingTime :: (HasLock e, HasLogHandles e, HasStartingTime e, HasVerbosity e) => Logger e ()
logInfoStartingTime = do
  ti <- reader getStartingTime
  logInfoS $ "Starting time: " <> renderTime ti

-- | Log end time.
logInfoEndTime :: (HasLock e, HasLogHandles e, HasStartingTime e, HasVerbosity e) => Logger e ()
logInfoEndTime = do
  ti <- reader getStartingTime
  te <- liftIO getCurrentTime
  let dt = te `diffUTCTime` ti
  logInfoB $ "Wall clock run time: " <> renderDuration dt <> "."
  logInfoS $ "End time: " <> renderTime te
