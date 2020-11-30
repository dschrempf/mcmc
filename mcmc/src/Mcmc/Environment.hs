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
import System.IO

-- | The environment of an MCMC run.
data Environment = Environment
  { settings :: Settings,
    -- | We have to use 'Maybe' here, because we do not want to open any log
    -- file when being 'Quiet'.
    logHandle :: Maybe Handle,
    -- | Number of used capabilities.
    numCap :: Int,
    -- | Used to calculate the ETA.
    startingTime :: UTCTime
  }
  deriving (Eq, Show)

-- | Initialize the environment.
--
-- Open log file, get current time.
initializeEnvironment ::
  Settings ->
  -- | Number of capabilities.
  Int ->
  IO Environment
initializeEnvironment s c = do
  t <- getCurrentTime
  mh <- case sVerbosity s of
    Quiet -> return Nothing
    _ -> do
      h <- openWithExecutionMode em fn
      return $ Just h
  return $ Environment s mh c t
  where
    fn = fromAnalysisName (sAnalysisName s) ++ ".log"
    em = sExecutionMode s
