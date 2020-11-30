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
    aParallelizationCheck,
  )
where

import Control.Concurrent
import qualified Data.ByteString.Lazy.Char8 as BL
import Mcmc.Environment

-- | Class for algorithms used by MCMC samplers.
class Algorithm a where
  aName :: a -> String

  aIteration :: a -> Int

  -- TODO: Splitmix. Remove IO monad constraint as soon as possible.

  aIterate ::
    -- | Number of capabilities for parallel execution.
    Int ->
    a ->
    IO a

  aAutoTune :: a -> a

  aResetAcceptance :: a -> a

  aSummarizeCycle :: a -> BL.ByteString

  -- | Open all monitor files and provide the file handles.
  aOpenMonitors :: Environment -> a -> IO a

  -- TODO: Can I remove the 'Environment' argument? See 'mc3ExecuteMonitors'.

  -- | Execute file monitors and possible return a string to be written to the
  -- standard output and the log file.
  aExecuteMonitors :: Environment -> a -> IO (Maybe BL.ByteString)

  -- | Header of monitor to standard output.
  aStdMonitorHeader :: a -> BL.ByteString

  -- | Close all files and remove the file handles.
  aCloseMonitors :: a -> IO a

  -- | Save chain(s).
  aSave ::
    -- | Maximum length of trace.
    Int ->
    -- | Analysis name.
    String ->
    a ->
    IO ()

-- TODO: Splitmix. Guess what? Remove IO.
-- | For a given algorithm check if parallelization is beneficial.
aParallelizationCheck ::
  Algorithm a =>
  a ->
  -- Number of capabilities to use with 'aIterate'.
  IO Int
aParallelizationCheck a = do
  c <- getNumCapabilities
  -- TODO.
  undefined
