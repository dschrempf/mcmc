{-# LANGUAGE MultiParamTypeClasses #-}

-- |
-- Module      :  Mcmc.Algorithm
-- Description :  Algortihms for Markov chain Monte Carlo samplers
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 16 14:37:11 2020.
module Mcmc.Algorithm
  ( Algorithm (..),
  )
where

import qualified Data.ByteString.Lazy.Char8 as BL
import Mcmc.Environment
import Numeric.Log

-- | @t@ is the algorithm, @a@ is the state space.
class Algorithm t a where
  aName :: t a -> String

  aIteration :: t a -> Int

  -- TODO: Splitmix. Remove IO monad as soon as possible.

  aIterate :: t a -> IO (t a)

  aAutoTune :: t a -> t a

  aResetAcceptance :: t a -> t a

  aSummarizeCycle :: t a -> BL.ByteString

  -- | Open all monitor files and provide the file handles.
  aOpenMonitors :: Environment -> t a -> IO (t a)

  -- | Execute file monitors and possible return a string to be written to the
  -- standard output and the log file.
  aExecuteMonitors :: Environment -> t a -> IO (Maybe BL.ByteString)

  -- | Header of monitor to standard output.
  aStdMonitorHeader :: t a -> BL.ByteString

  -- | Close all files and remove the file handles.
  aCloseMonitors :: t a -> IO (t a)

  -- | Save chain(s) with trace of given maximum length.
  aSave ::
    Int ->
    -- | Analysis name.
    String ->
    t a ->
    IO ()
