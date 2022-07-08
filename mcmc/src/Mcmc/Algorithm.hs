-- |
-- Module      :  Mcmc.Algorithm
-- Description :  Algortihms for Markov chain Monte Carlo samplers
-- Copyright   :  2021 Dominik Schrempf
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
import Mcmc.Cycle
import Mcmc.Proposal
import Mcmc.Settings

-- | Class for algorithms used by MCMC samplers.
class Algorithm a where
  -- | Name.
  aName :: a -> String

  -- | Current iteration.
  aIteration :: a -> Int

  -- | Check if the current state is invalid. At the moment this should just
  -- check if the posterior probability is zero or NaN.
  aIsInValidState :: a -> Bool

  -- | Sample the next state.
  aIterate :: IterationMode -> ParallelizationMode -> a -> IO a

  -- | Auto tune all proposals over the last N iterations.
  --
  -- NOTE: Computation in the 'IO' Monad is necessary because the trace is
  -- mutable.
  aAutoTune :: TuningType -> Int -> a -> IO a

  -- | Reset acceptance counts.
  aResetAcceptance :: a -> a

  -- | Clean after burn in. In particular, this is used to reduce the length of
  -- the trace, if required.
  aCleanAfterBurnIn :: TraceLength -> a -> IO a

  -- | Summarize the cycle.
  aSummarizeCycle :: IterationMode -> a -> BL.ByteString

  -- | Open required monitor files and setup corresponding file handles.
  aOpenMonitors :: AnalysisName -> ExecutionMode -> a -> IO a

  -- | Execute file monitors and possibly return a monitor string to be written
  -- to the standard output and the log file.
  aExecuteMonitors ::
    Verbosity ->
    -- | Starting time.
    UTCTime ->
    -- | Total number of iterations including burn in.
    Int ->
    a ->
    IO (Maybe BL.ByteString)

  -- | Header of monitor to standard output.
  aStdMonitorHeader :: a -> BL.ByteString

  -- | Close monitor files and remove the file handles.
  aCloseMonitors :: a -> IO a

  -- | Save analysis.
  aSave :: AnalysisName -> a -> IO ()
