{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Settings
-- Description :  Settings of the Markov chain Monte Carlo Sampler
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 16 11:13:01 2020.
module Mcmc.Settings
  ( BurnIn (..),
    ExecutionMode (..),
    openWithExecutionMode,
    SaveMode (..),
    Verbosity (..),
    Settings (..),
    loadSettings,
  )
where

import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Lazy.Char8 as BL
import System.Directory
import System.IO

-- | Burn in specification.
data BurnIn
  = -- | No burn in.
    NoBurnIn
  | -- | Burn in for a given number of iterations.
    BurnInNoAutoTuning Int
  | -- | Burn in for a given number of iterations. Auto tuning with a given auto
    -- tuning period is enabled.
    BurnInWithAutoTuning Int Int
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''BurnIn)

-- | Execution mode.
data ExecutionMode
  = -- | Call 'error' if an output files exists.
    Fail
  | -- | Overwrite existing output files.
    Overwrite
  | -- | Continue a previous run and append to output files. Call 'error' if an
    -- output file does not exist.
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

-- | Settings of the Markov chain Monte Carlo sampler.
data Settings = Settings
  { -- | Name of the Markov chain Monte Carlo sampler.
    name :: String,
    burnIn :: BurnIn,
    -- | Number of normal iterations excluding burn in. Note that auto tuning
    -- only happens during burn in.
    iterations :: Int,
    executionMode :: ExecutionMode,
    saveMode :: SaveMode,
    verbosity :: Verbosity
  }
  deriving (Eq, Show)

$(deriveJSON defaultOptions ''Settings)

-- | Load settings from a file.
loadSettings :: FilePath -> IO Settings
loadSettings fn = either error id . eitherDecode <$> BL.readFile fn
