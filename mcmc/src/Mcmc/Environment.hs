-- |
-- Module      :  Mcmc.Environment
-- Description :  Runtime environment
-- Copyright   :  2021 Dominik Schrempf
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
    closeEnvironment,
  )
where

import Control.Concurrent.MVar
import Control.Monad
import Data.Time
import Mcmc.Logger
import Mcmc.Settings
import System.IO

-- | The environment of an MCMC run.
data Environment s = Environment
  { settings :: s,
    -- | List will be empty if using 'Quiet'. If 'LogStdOutAndFile' is used
    -- 'logHandles' contains two handles to the standard output and the log
    -- file.
    logHandles :: [Handle],
    -- | MVar blocking output.
    outLock :: MVar (),
    -- | Used to calculate the ETA.
    startingTime :: UTCTime
  }
  deriving (Eq)

instance HasExecutionMode s => HasExecutionMode (Environment s) where
  getExecutionMode = getExecutionMode . settings

instance HasLock (Environment s) where
  getLock = outLock

instance HasLogHandles (Environment s) where
  getLogHandles = logHandles

instance HasStartingTime (Environment s) where
  getStartingTime = startingTime

instance HasLogMode s => HasLogMode (Environment s) where
  getLogMode = getLogMode . settings

instance HasVerbosity s => HasVerbosity (Environment s) where
  getVerbosity = getVerbosity . settings

-- | Initialize the environment.
--
-- Open log file, get current time.
initializeEnvironment ::
  (HasAnalysisName s, HasExecutionMode s, HasLogMode s, HasVerbosity s) =>
  s ->
  IO (Environment s)
initializeEnvironment s = do
  t <- getCurrentTime
  mh <- case (getLogMode s, getVerbosity s) of
    (_, Quiet) -> return []
    (LogStdOutAndFile, _) -> do
      h <- openWithExecutionMode em fn
      return [stdout, h]
    (LogFileOnly, _) -> do
      h <- openWithExecutionMode em fn
      return [h]
    (LogStdOutOnly, _) -> return [stdout]
  lock <- newMVar ()
  return $ Environment s mh lock t
  where
    fn = fromAnalysisName (getAnalysisName s) ++ ".mcmc.log"
    em = getExecutionMode s

-- | Close file handles.
closeEnvironment :: Environment s -> IO ()
closeEnvironment e = forM_ hs hClose
  where
    hs = filter (/= stdout) $ logHandles e
