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
import Mcmc.Logger
import Mcmc.Settings
import System.IO

-- | The environment of an MCMC run.
data Environment s = Environment
  { settings :: s,
    -- | We have to use 'Maybe' here, because we do not want to open any log
    -- file when being 'Quiet'.
    logHandle :: Maybe Handle,
    -- | Used to calculate the ETA.
    startingTime :: UTCTime
  }
  deriving (Eq, Show)

instance HasExecutionMode s => HasExecutionMode (Environment s) where
  getExecutionMode = getExecutionMode . settings

instance HasMaybeLogHandle (Environment s) where
  getMaybeLogHandle = logHandle

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
  return $ Environment s mh t
  where
    fn = fromAnalysisName (getAnalysisName s) ++ ".mcmc.log"
    em = getExecutionMode s
