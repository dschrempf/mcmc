{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Settings
-- Description :  Settings of Markov chain Monte Carlo samplers
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 16 11:13:01 2020.
module Mcmc.Settings
  ( -- * Data types
    AnalysisName (..),
    BurnInSpecification (..),
    burnInIterations,
    Iterations (..),
    ExecutionMode (..),
    openWithExecutionMode,
    ParallelizationMode (..),
    SaveMode (..),
    Verbosity (..),

    -- * Settings
    Settings (..),
    settingsSave,
    settingsLoad,
    settingsCheck,
  )
where

import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Lazy.Char8 as BL
import System.Directory
import System.IO

-- | Analysis name of the MCMC sampler.
newtype AnalysisName = AnalysisName {fromAnalysisName :: String}
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''AnalysisName)

-- | Burn in specification.
data BurnInSpecification
  = -- | No burn in.
    NoBurnIn
  | -- | Burn in for a given number of iterations.
    BurnInWithoutAutoTuning Int
  | -- | Burn in for a given number of iterations. Enable auto tuning with a
    -- given period.
    BurnInWithAutoTuning Int Int
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''BurnInSpecification)

-- | Get the number of burn in iterations.
burnInIterations :: BurnInSpecification -> Int
burnInIterations NoBurnIn = 0
burnInIterations (BurnInWithoutAutoTuning n) = n
burnInIterations (BurnInWithAutoTuning n _) = n

-- | Number of normal iterations after burn in.
--
-- Note that auto tuning only happens during burn in.
newtype Iterations = Iterations {fromIterations :: Int}
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''Iterations)

-- | Execution mode.
data ExecutionMode
  = -- | Perform new run.
    --
    -- Call 'error' if an output files exists.
    Fail
  | -- | Perform new run.
    --
    -- Overwrite existing output files.
    Overwrite
  | -- | Continue a previous run and append to output files.
    --
    -- Call 'error' if an output file does not exist.
    Continue
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''ExecutionMode)

-- | Open a file honoring the execution mode.
--
-- Call 'error' if execution mode is
--
-- - 'Continue' and file does not exist.
--
-- - 'Fail' and file exists.
openWithExecutionMode :: ExecutionMode -> FilePath -> IO Handle
openWithExecutionMode em fn = do
  fe <- doesFileExist fn
  case (em, fe) of
    (Continue, False) ->
      error $ "openWithExecutionMode: Cannot continue; file does not exist: " ++ fn ++ "."
    (Continue, True) ->
      openFile fn AppendMode
    (Fail, True) ->
      error $ "openWithExecutionMode: File exists: " ++ fn ++ "; use 'Overwrite'?"
    _ -> openFile fn WriteMode

-- | Parallelization mode.
--
-- Parallel execution of the chains is only beneficial when the algorithm allows
-- for parallelization, and if computation of the next iteration takes a long
-- time. If the calculation of the next state is fast, sequential execution is
-- usually beneficial, even for algorithms involving parallel chains. If the
-- calculation of the next state is slow, parallel execution may be beneficial.
--
-- - The "Mcmc.Algorithm.Metropolis" algorithm is inherently sequential.
--
-- - The "Mcmc.Algorithm.MC3" algorithm works well with parallelization.
--
-- Of course, also the prior or likelihood functions can be computed in
-- parallel. However, this library is not aware of how these functions are
-- computed.
data ParallelizationMode
  = Sequential
  | Parallel
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''ParallelizationMode)

-- | Should the MCMC run with trace of given maximum length be saved at the end
-- of the run?
data SaveMode = NoSave | SaveWithTrace Int
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SaveMode)

-- $(deriveJSON defaultOptions ''CleaningMode)

-- | Not much to say here.
data Verbosity = Quiet | Warn | Info | Debug
  deriving (Eq, Ord, Read, Show)

$(deriveJSON defaultOptions ''Verbosity)

-- | Settings of an MCMC sampler.
data Settings = Settings
  { sAnalysisName :: AnalysisName,
    sBurnIn :: BurnInSpecification,
    sIterations :: Iterations,
    sExecutionMode :: ExecutionMode,
    sParallelizationMode :: ParallelizationMode,
    sSaveMode :: SaveMode,
    sVerbosity :: Verbosity
  }
  deriving (Eq, Show)

$(deriveJSON defaultOptions ''Settings)

settingsFn :: String -> FilePath
settingsFn n = n ++ ".settings"

-- | Save settings to a file determined by the analysis name.
settingsSave :: Settings -> IO ()
settingsSave s = BL.writeFile fn $ encode s
  where
    fn = settingsFn $ fromAnalysisName $ sAnalysisName s

-- | Load settings.
settingsLoad :: AnalysisName -> IO Settings
settingsLoad (AnalysisName n) = either error id . eitherDecode <$> BL.readFile fn
  where
    fn = settingsFn n

-- Show settings and call 'error'.
settingsError :: Settings -> Int -> String -> a
settingsError s i err =
  error $
    show s
      ++ "\n"
      ++ "Current iteration: "
      ++ show i
      ++ "\n"
      ++ "settingsError: "
      ++ err

-- | Check settings.
--
-- Call 'error' if:
--
-- - The analysis name is the empty string.
--
-- - The number of burn in iterations is negative.
--
-- - Auto tuning period is zero or negative.
--
-- - The number of iterations is negative.
--
-- - The current iteration is larger than the total number of iterations.
--
-- - The current iteration is non-zero but the execution mode is not 'Continue'.
--
-- - The current iteration is zero but the execution mode is 'Continue'.
settingsCheck ::
  Settings ->
  -- | Current iteration.
  Int ->
  IO ()
settingsCheck s@(Settings nm bi i em _ _ _) iCurrent
  | null (fromAnalysisName nm) = serr "Analysis name is the empty string."
  | burnInIterations bi < 0 = serr "Number of burn in iterations is negative."
  | not $ burnInAutoTuningPeriodValid bi = serr "Auto tuning period is zero or negative."
  | fromIterations i < 0 = serr "Number of iterations is negative."
  | burnInIterations bi + fromIterations i - iCurrent < 0 =
    serr "Current iteration is larger than the total number of iterations."
  | iCurrent /= 0 && em /= Continue =
    serr "Current iteration is non-zero but execution mode is not 'Continue'."
  | iCurrent == 0 && em == Continue =
    serr "Current iteration is zero but execution mode is 'Continue'."
  | otherwise = return ()
  where
    serr = settingsError s iCurrent
    burnInAutoTuningPeriodValid :: BurnInSpecification -> Bool
    burnInAutoTuningPeriodValid (BurnInWithAutoTuning _ t) = t > 0
    burnInAutoTuningPeriodValid _ = True
