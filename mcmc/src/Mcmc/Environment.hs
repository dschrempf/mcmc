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

import Control.Concurrent.MVar
import Data.Time.Clock
import Mcmc.Logger
import Mcmc.Settings
import System.IO

-- | The environment of an MCMC run.
data Environment s = Environment
  { settings :: s,
    -- | We have to use 'Maybe' here, because we do not want to open any log
    -- file when being 'Quiet'.
    logHandle :: Maybe Handle,
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

instance HasMaybeLogHandle (Environment s) where
  getMaybeLogHandle = logHandle

instance HasStartingTime (Environment s) where
  getStartingTime = startingTime

instance HasVerbosity s => HasVerbosity (Environment s) where
  getVerbosity = getVerbosity . settings

-- | Initialize the environment.
--
-- Open log file, get current time.
initializeEnvironment ::
  (HasAnalysisName s, HasExecutionMode s, HasVerbosity s) =>
  s ->
  IO (Environment s)
initializeEnvironment s = do
  t <- getCurrentTime
  mh <- case getVerbosity s of
    Quiet -> return Nothing
    _ -> do
      h <- openWithExecutionMode em fn
      return $ Just h
  lock <- newMVar ()
  return $ Environment s mh lock t
  where
    fn = fromAnalysisName (getAnalysisName s) ++ ".mcmc.log"
    em = getExecutionMode s
