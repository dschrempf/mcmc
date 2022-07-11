{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Logger
-- Description :  Minimal monad logger
-- Copyright   :  2021 Dominik Schrempf
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
import Control.Monad.Reader
import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Time.Clock
import Data.Version (showVersion)
import Mcmc.Monitor.Time
import Paths_mcmc (version)
import System.IO

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

msgPrepare :: BL.ByteString -> BL.ByteString -> BL.ByteString
msgPrepare pref msg = BL.intercalate "\n" $ map (BL.append pref) $ BL.lines msg

-- Make sure that concurrent output is not scrambled.
atomicAction :: (HasLock e, MonadReader e m, MonadIO m) => IO () -> m ()
atomicAction a = do
  l <- asks getLock
  liftIO $ withMVar l (const a)

-- | Write to standard output and maybe to log file.
logOutB ::
  (HasLogHandles e, HasLock e, MonadReader e m, MonadIO m) =>
  -- | Prefix.
  BL.ByteString ->
  -- | Message.
  BL.ByteString ->
  m ()
logOutB pref msg = do
  hs <- asks getLogHandles
  mapM_ (atomicAction . (`BL.hPutStrLn` msg')) hs
  where
    msg' = msgPrepare pref msg

-- Perform debug action.
logDebugA :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => m () -> m ()
logDebugA a = asks getVerbosity >>= \v -> when (v >= Debug) a

-- | Log debug message.
logDebugB :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => BL.ByteString -> m ()
logDebugB = logDebugA . logOutB "D: "

-- | Log debug message.
logDebugS :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => String -> m ()
logDebugS = logDebugB . BL.pack

-- Perform warning action.
logWarnA :: (HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => m () -> m ()
logWarnA a = asks getVerbosity >>= \v -> when (v >= Warn) a

-- | Log warning message.
logWarnB :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => BL.ByteString -> m ()
logWarnB = logWarnA . logOutB "W: "

-- | Log warning message.
logWarnS :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => String -> m ()
logWarnS = logWarnB . BL.pack

-- Perform info action.
logInfoA :: (HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => m () -> m ()
logInfoA a = asks getVerbosity >>= \v -> when (v >= Info) a

-- | Log info message.
logInfoB :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => BL.ByteString -> m ()
logInfoB = logInfoA . logOutB "   "

-- | Log info message.
logInfoS :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => String -> m ()
logInfoS = logInfoB . BL.pack

-- | Log info header.
logInfoHeader :: (HasLock e, HasLogHandles e, HasVerbosity e, MonadReader e m, MonadIO m) => m ()
logInfoHeader = do
  logInfoS (replicate 70 '-')
  logInfoS ("MCMC sampler; version " ++ showVersion version <> ".")
  logInfoS "Developed by: Dominik Schrempf."
  logInfoS "License: GPL-3.0-or-later."
  logInfoS (replicate 70 '-')

-- | Log starting time.
logInfoStartingTime ::
  ( HasLock e,
    HasLogHandles e,
    HasStartingTime e,
    HasVerbosity e,
    MonadReader e m,
    MonadIO m
  ) =>
  m ()
logInfoStartingTime = do
  ti <- asks getStartingTime
  logInfoS $ "Starting time: " <> renderTime ti

-- | Log end time.
logInfoEndTime ::
  ( HasLock e,
    HasLogHandles e,
    HasStartingTime e,
    HasVerbosity e,
    MonadReader e m,
    MonadIO m
  ) =>
  m ()
logInfoEndTime = do
  ti <- asks getStartingTime
  te <- liftIO getCurrentTime
  let dt = te `diffUTCTime` ti
  logInfoB $ "Wall clock run time: " <> renderDuration dt <> "."
  logInfoS $ "End time: " <> renderTime te
