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
  ( Overwrite (..),
    SaveChain (..),
    Verbosity (..),
    Environment (..),
    forceOverwrite,
    saveN,
    quiet,
    debug,
  )
where

import Data.Aeson.TH
import Data.Default

-- | Force overwrite of output files, or fail with an error message?
data Overwrite = Fail | Force
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''Overwrite)

-- | Save the chain with trace of given maximum length at the end of the run?
data SaveChain = NoSave | SaveN Int
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SaveChain)

-- | Not much to say here.
data Verbosity = Quiet | Warn | Info | Debug
  deriving (Eq, Ord, Read, Show)

$(deriveJSON defaultOptions ''Verbosity)

-- | Environment of the Markov chain Monte Carlo sampler.
data Environment = Environment
  { overwrite :: Overwrite,
    saveChain :: SaveChain,
    verbosity :: Verbosity
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''Environment)

instance Default Environment where def = Environment Fail NoSave Info

-- | Force overwrite of output files.
forceOverwrite :: Environment -> Environment
forceOverwrite e = e {overwrite = Force}

-- | Save the chain with trace of given maximum length at the end of the run.
saveN :: Int -> Environment -> Environment
saveN n e = e {saveChain = SaveN n}

-- | Be quiet.
quiet :: Environment -> Environment
quiet e = e {verbosity = Quiet}

-- | Show debug output.
debug :: Environment -> Environment
debug e = e {verbosity = Debug}
