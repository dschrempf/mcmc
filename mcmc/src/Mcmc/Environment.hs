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
    logHandle :: Maybe Handle,
    -- | Used to calculate the ETA.
    startingTime :: UTCTime
  }
  deriving (Eq, Show)

-- | Initialize the environment.
--
-- Open log file, get current time.
initializeEnvironment :: Settings -> IO Environment
initializeEnvironment s = do
  t <- getCurrentTime
  case sVerbosity s of
    Quiet -> return $ Environment s Nothing t
    _ -> do
      h <- openWithExecutionMode em fn
      return $ Environment s (Just h) t
  where
    fn = sAnalysisName s ++ ".log"
    em = sExecutionMode s
