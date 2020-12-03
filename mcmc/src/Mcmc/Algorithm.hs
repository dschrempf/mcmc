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
import Data.Time
import Mcmc.Settings

-- | Class for algorithms used by MCMC samplers.
class Algorithm a where
  aName :: a -> String

  aIteration :: a -> Int

  -- TODO: Splitmix. Remove IO monad constraint as soon as possible.

  aIterate ::
    ParallelizationMode ->
    a ->
    IO a

  aAutoTune :: a -> a

  aResetAcceptance :: a -> a

  aSummarizeCycle :: a -> BL.ByteString

  -- | Open all monitor files and provide the file handles.
  aOpenMonitors :: AnalysisName -> ExecutionMode -> a -> IO a

  -- | Execute file monitors and possible return a string to be written to the
  -- standard output and the log file.
  aExecuteMonitors ::
    Verbosity ->
    -- | Starting time.
    UTCTime ->
    -- | Total number of iterations.
    Int ->
    a ->
    IO (Maybe BL.ByteString)

  -- | Header of monitor to standard output.
  aStdMonitorHeader :: a -> BL.ByteString

  -- | Close all files and remove the file handles.
  aCloseMonitors :: a -> IO a

  -- TODO: Splitmix. Guess what? Remove IO. But then, we remain with PrimMonad
  -- because of the trace.

  -- | Save chain(s).
  aSave ::
    AnalysisName ->
    a ->
    IO ()

