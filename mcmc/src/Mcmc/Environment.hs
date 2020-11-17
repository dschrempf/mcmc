-- |
-- Module      :  Mcmc.Environment
-- Description :  Runtime environment
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Nov 17 11:00:09 2020.
module Mcmc.Environment
  ( Environment (..),
    initializeEnvironment,
  )
where

import Data.Time.Clock
import Mcmc.Settings
import System.Directory
import System.IO

data Environment = Environment
  { settings :: Settings,
    -- | The log handle is set internally with 'openLogFile'.
    logHandle :: Handle,
    -- | Starting time. Used to calculate the ETA.
    startingTime :: UTCTime
  }
  deriving (Eq, Show)

-- | Initialize the environment.
--
-- Open log file, get current time.
--
-- Call 'error' if the log file exists and output mode is 'Fail'.
initializeEnvironment :: Settings -> IO Environment
initializeEnvironment s = do
  fe <- doesFileExist fn
  h <- case (fe, om) of
    (False, _) -> openFile fn WriteMode
    (True, Overwrite) -> openFile fn WriteMode
    (True, Fail) -> error "openLogFile: Log file exists."
    (True, Append) -> openFile fn AppendMode
  t <- getCurrentTime
  return $ Environment s h t
  where
    nm = name s
    fn = nm ++ ".log"
    om = outputMode s
