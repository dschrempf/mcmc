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
module Mcmc.Internal.Logger
  ( Verbosity (..),
    LogEnv (..),
    logOutB,
    logDebugA,
    logDebugB,
    logDebugS,
    logWarnA,
    logWarnB,
    logWarnS,
    logInfoA,
    logInfoB,
    logInfoS,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Reader
import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Lazy.Char8 as BL
import System.IO

-- | Not much to say here.
data Verbosity = Quiet | Warn | Info | Debug
  deriving (Eq, Ord, Read, Show)

$(deriveJSON defaultOptions ''Verbosity)

-- | Environment with logging information.
class LogEnv e where
  getMaybeLogHandle :: e -> Maybe Handle
  getVerbosity :: e -> Verbosity

-- | Reader transformer used for logging to a file and to standard output.
type Logger e a = ReaderT e IO a

msgPrepare :: BL.ByteString -> BL.ByteString -> BL.ByteString
msgPrepare pref msg = BL.intercalate "\n" $ map (BL.append pref) $ BL.lines msg

-- | Write to standard output and maybe to log file.
logOutB ::
  LogEnv e =>
  -- | Prefix.
  BL.ByteString ->
  -- | Message.
  BL.ByteString ->
  Logger e ()
logOutB pref msg = do
  liftIO $ BL.putStrLn msg'
  mh <- reader getMaybeLogHandle
  case mh of
    Nothing -> return ()
    Just h -> liftIO $ BL.hPutStrLn h msg'
  where
    msg' = msgPrepare pref msg

-- | Perform debug action.
logDebugA :: LogEnv e => Logger e () -> Logger e ()
logDebugA a = reader getVerbosity >>= \v -> when (v == Debug) a

-- | Print debug message.
logDebugB :: LogEnv e => BL.ByteString -> Logger e ()
logDebugB = logDebugA . logOutB "D: "

-- | Print debug message.
logDebugS :: LogEnv e => String -> Logger e ()
logDebugS = logDebugB . BL.pack

-- | Perform warning action.
logWarnA :: LogEnv e => Logger e () -> Logger e ()
logWarnA a = reader getVerbosity >>= \v -> when (v >= Warn) a

-- | Print warning message.
logWarnB :: LogEnv e => BL.ByteString -> Logger e ()
logWarnB = logWarnA . logOutB "W: "

-- | Print warning message.
logWarnS :: LogEnv e => String -> Logger e ()
logWarnS = logWarnB . BL.pack

-- | Perform info action.
logInfoA :: LogEnv e => Logger e () -> Logger e ()
logInfoA a = reader getVerbosity >>= \v -> when (v >= Info) a

-- | Print info message.
logInfoB :: LogEnv e => BL.ByteString -> Logger e ()
logInfoB = logInfoA . logOutB "I: "

-- | Print info message.
logInfoS :: LogEnv e => String -> Logger e ()
logInfoS = logInfoB . BL.pack
