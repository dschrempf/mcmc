{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Verbosity
-- Description :  Be quiet! Or better not?
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Sat Jun 27 10:49:28 2020.
module Mcmc.Verbosity
  ( Verbosity (..),
    warn,
    info,
    debug,
  )
where

import Control.Monad
import Data.Aeson.TH

-- | Not much to say here.
data Verbosity = Quiet | Warn | Info | Debug deriving (Show, Eq, Ord)

$(deriveJSON defaultOptions ''Verbosity)

-- | Perform action if 'Verbosity' is 'Warn' or higher.
warn :: Applicative m => Verbosity -> m () -> m ()
warn v = when (v >= Warn)

-- | Perform action if 'Verbosity' is 'Info' or higher.
info :: Applicative m => Verbosity -> m () -> m ()
info v = when (v >= Info)

-- | Perform action if 'Verbosity' is 'Debug'.
debug :: Applicative m => Verbosity -> m () -> m ()
debug v = when (v == Debug)
