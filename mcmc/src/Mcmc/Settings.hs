{-# LANGUAGE DerivingVia #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Settings
-- Description :  Settings of Markov chain Monte Carlo samplers
-- Copyright   :  (c) Dominik Schrempf, 2021
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
    HasAnalysisName (..),
    BurnInSpecification (..),
    burnInIterations,
    Iterations (..),
    TraceLength (..),
    ExecutionMode (..),
    HasExecutionMode (..),
    openWithExecutionMode,
    ParallelizationMode (..),
    SaveMode (..),
    LogMode (..),
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
import Mcmc.Logger
import System.Directory
import System.IO

-- | Analysis name of the MCMC sampler.
newtype AnalysisName = AnalysisName {fromAnalysisName :: String}
  deriving (Eq, Read, Show)
  deriving (Monoid, Semigroup) via String

$(deriveJSON defaultOptions ''AnalysisName)

-- | Types with analysis names.
class HasAnalysisName s where
  getAnalysisName :: s -> AnalysisName

-- | Burn in specification.
data BurnInSpecification
  = -- | No burn in.
    NoBurnIn
  | -- | Burn in for a given number of iterations.
    BurnInWithoutAutoTuning Int
  | -- | Burn in for a given number of iterations. Enable auto tuning with a
    -- given period.
    BurnInWithAutoTuning Int Int
  | -- | Burn in with the given list of auto tuning periods.
    --
    -- For example, @BurnInWithCustomAutoTuning [100,200]@ performs 300
    -- iterations with two auto tuning steps. One after 100 iterations, the
    -- second one after 200 more iterations.
    --
    -- Usually it is useful to auto tune more frequently in the beginning of the
    -- MCMC run.
    BurnInWithCustomAutoTuning [Int]
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''BurnInSpecification)

-- Check if the burn in specification is valid.
burnInValid :: BurnInSpecification -> Bool
burnInValid NoBurnIn = True
burnInValid (BurnInWithoutAutoTuning n) = n > 0
burnInValid (BurnInWithAutoTuning n t) = n > 0 && t > 0
burnInValid (BurnInWithCustomAutoTuning xs) = not (null xs) && all (> 0) xs

-- | Get the number of burn in iterations.
burnInIterations :: BurnInSpecification -> Int
burnInIterations NoBurnIn = 0
burnInIterations (BurnInWithoutAutoTuning n) = n
burnInIterations (BurnInWithAutoTuning n _) = n
burnInIterations (BurnInWithCustomAutoTuning xs) = sum xs

-- | Number of normal iterations after burn in.
--
-- Note that auto tuning only happens during burn in.
newtype Iterations = Iterations {fromIterations :: Int}
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''Iterations)

-- | The length of the stored "Mcmc.Chain.Trace".
--
-- Be careful, this setting determines the memory requirement of the MCMC chain.
data TraceLength
  = -- | Automatically determine the length of the trace. The value is
    -- determined by the 'Mcmc.Monitor.MonitorBatch' with largest batch size.
    TraceAuto
  | -- | Store a given minimum number of iterations of the chain. Store more
    --  iterations if required (see 'TraceAuto').
    TraceMinimum Int

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

-- | Types with execution modes.
class HasExecutionMode s where
  getExecutionMode :: s -> ExecutionMode

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
    _ -> do
      h <- openFile fn WriteMode
      hSetBuffering h LineBuffering
      return h

-- One could automatically select 'Parallel' or 'Sequential' according to the
-- number of capabilities when initializing the environment or according to the
-- iteration time in dependence of the number of used capabilities. However, I
-- decided to opt for a manual configuration, because more capabilities may be
-- available and other parts of the program may be executed in parallel even if
-- sequential execution of the MCMC sampler is beneficial.

-- | Parallelization mode.
--
-- Parallel execution of the chains is only beneficial when the algorithm allows
-- for parallelization, and if computation of the next iteration takes some
-- time. If the calculation of the next state is fast, sequential execution is
-- usually beneficial, even for algorithms involving parallel chains.
--
-- - The "Mcmc.Algorithm.MHG" algorithm is inherently sequential.
--
-- - The "Mcmc.Algorithm.MC3" algorithm works well with parallelization.
--
-- Of course, also the prior or likelihood functions can be computed in
-- parallel. However, this library is unaware about how these functions are
-- computed.
data ParallelizationMode
  = Sequential
  | Parallel
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''ParallelizationMode)

-- | Define information stored on disk.
data SaveMode
  = -- | Do not save the MCMC analysis. The analysis can not be continued.
    NoSave
  | -- | Save the MCMC analysis so that it can be continued. This can be slow,
    -- if the trace is long, or if the states are large objects. See
    -- 'TraceLength'.
    Save
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SaveMode)

-- | Settings of an MCMC sampler.
data Settings = Settings
  { sAnalysisName :: AnalysisName,
    sBurnIn :: BurnInSpecification,
    sIterations :: Iterations,
    sExecutionMode :: ExecutionMode,
    sParallelizationMode :: ParallelizationMode,
    sSaveMode :: SaveMode,
    sLogMode :: LogMode,
    sVerbosity :: Verbosity
  }
  deriving (Eq, Show)

instance HasAnalysisName Settings where
  getAnalysisName = sAnalysisName

instance HasExecutionMode Settings where
  getExecutionMode = sExecutionMode

instance HasLogMode Settings where
  getLogMode = sLogMode

instance HasVerbosity Settings where
  getVerbosity = sVerbosity

$(deriveJSON defaultOptions ''Settings)

settingsFn :: String -> FilePath
settingsFn n = n ++ ".mcmc.settings"

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
settingsCheck s@(Settings nm bi i em _ _ _ _) iCurrent
  | null (fromAnalysisName nm) = serr "Analysis name is the empty string."
  | burnInIterations bi < 0 = serr "Number of burn in iterations is negative."
  | not $ burnInValid bi = serr $ "Burn in setting invalid: " <> show bi
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
