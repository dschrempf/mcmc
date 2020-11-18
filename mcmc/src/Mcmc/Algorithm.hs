-- |
-- Module      :  Mcmc.Algorithm
-- Description :  MCMC algorithms
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

-- | TODO: REFACTOR. Documentation.
class Algorithm a where
  algorithmName :: a -> String

  algorithmIteration :: a -> Int

  -- TODO: Splitmix. Remove IO monad as soon as possible.

  algorithmIterate :: a -> IO a

  algorithmAutoTune :: a -> a

  algorithmResetAcceptance :: a -> a

  algorithmSummarizeCycle :: a -> BL.ByteString

  -- | Open all monitor files and provide the file handles.
  algorithmOpenMonitors :: Environment -> a -> IO a

  -- | Execute file monitors and possible return a string to be written to the
  -- standard output and the log file.
  algorithmExecuteMonitors :: Environment -> a -> IO (Maybe BL.ByteString)

  -- | Close all files and remove the file handles.
  algorithmCloseMonitors :: a -> IO a

  -- | Save chain(s) with trace of given maximum length.
  algorithmSaveWith :: Int -> FilePath -> a -> IO ()

  -- | Report prior and likelihood; useful for debugging.
  algorithmReport :: a -> (Log Double, Log Double)
