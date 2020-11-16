{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Environment
-- Description :  Environment of the Markov chain Monte Carlo Sampler
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 16 11:13:01 2020.
module Mcmc.Environment
  ( BurnIn (..),
    OutputMode (..),
    SaveMode (..),
    Verbosity (..),
    Environment (..),
    openLogFile,
  )
where

-- TODO: REFACTOR. Check documentation.

import Data.Aeson
import Data.Aeson.TH
import Data.Maybe
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

-- | Overwrite output files, fail with an error message, or append?
data OutputMode = Overwrite | Fail | Append
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''OutputMode)

-- | Should the MCMC run with trace of given maximum length be saved at the end
-- of the run?
data SaveMode = NoSave | SaveWithTrace Int
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SaveMode)

-- | Not much to say here.
data Verbosity = Quiet | Warn | Info | Debug
  deriving (Eq, Ord, Read, Show)

$(deriveJSON defaultOptions ''Verbosity)

-- | Environment of the Markov chain Monte Carlo sampler; created with
-- 'environment'.
data Environment = Environment
  { -- | Name of the Markov chain Monte Carlo sampler.
    name :: String,
    burnIn :: BurnIn,
    -- | Number of normal iterations excluding burn in. Note that auto tuning
    -- only happens during burn in.
    iterations :: Int,
    outputMode :: OutputMode,
    saveMode :: SaveMode,
    verbosity :: Verbosity,
    -- | The log handle is set internally with 'openLogFile'.
    logHandle :: Maybe Handle
  }
  deriving (Eq, Show)

-- | The 'Handle' is not stored.
instance ToJSON Environment where
  toJSON (Environment nm bi is om sm vb _) =
    object
      [ "name" .= nm,
        "burnIn" .= bi,
        "iterations" .= is,
        "outputMode" .= om,
        "saveMode" .= sm,
        "verbosity" .= vb
      ]
  toEncoding (Environment nm bi is om sm vb _) =
    pairs $
        "name" .= nm
        <> "burnIn" .= bi
        <> "iterations" .= is
        <> "outputMode" .= om
        <> "saveMode" .= sm
        <> "verbosity" .= vb

-- | The 'Handle' is not restored.
instance FromJSON Environment where
  parseJSON = withObject "Environment" $ \v ->
    Environment
      <$> v .: "name"
      <*> v .: "burnIn"
      <*> v .: "iterations"
      <*> v .: "outputMode"
      <*> v .: "saveMode"
      <*> v .: "verbosity"
      <*> pure Nothing

-- | Open log file.
--
-- Call 'error' if:
--
-- - The log file has already been opened.
-- - The log file exists and output mode is 'Fail'.
openLogFile :: Environment -> IO Environment
openLogFile env
  | isJust $ logHandle env = error "openLogFile: Log file has already been opened"
  | otherwise = do
    fe <- doesFileExist fn
    h <- case (fe, om) of
      (False, _) -> openFile fn WriteMode
      (True, Overwrite) -> openFile fn WriteMode
      (True, Fail) -> error "openLogFile: Log file exists."
      (True, Append) -> openFile fn AppendMode
    return $ env {logHandle = Just h}
  where
    nm = name env
    fn = nm ++ ".log"
    om = outputMode env

-- TODO.
-- mcmcDebugS $ "Log file name: " ++ lfn ++ "."
-- mcmcDebugB "Log file opened."
