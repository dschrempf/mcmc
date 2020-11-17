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
    OutputMode (..),
    SaveMode (..),
    -- CleaningMode (..),
    Verbosity (..),
    Settings (..),
  )
where

-- TODO: REFACTOR. Check documentation.

import Data.Aeson
import Data.Aeson.TH

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

-- -- | Should the state be cleaned periodically? A 'Mcmc.Chain.Cleaner' has to be provided
-- data CleaningMode = NoClean | CleanEvery Int
--   deriving (Eq, Read, Show)

-- $(deriveJSON defaultOptions ''CleaningMode)

-- | Not much to say here.
data Verbosity = Quiet | Warn | Info | Debug
  deriving (Eq, Ord, Read, Show)

$(deriveJSON defaultOptions ''Verbosity)

-- | Settings of the Markov chain Monte Carlo sampler; created with
-- 'environment'.
data Settings = Settings
  { -- | Name of the Markov chain Monte Carlo sampler.
    name :: String,
    burnIn :: BurnIn,
    -- | Number of normal iterations excluding burn in. Note that auto tuning
    -- only happens during burn in.
    iterations :: Int,
    outputMode :: OutputMode,
    saveMode :: SaveMode,
    -- cleaningMode :: CleaningMode,
    verbosity :: Verbosity
  }
  deriving (Eq, Show)

$(deriveJSON defaultOptions ''Settings)
