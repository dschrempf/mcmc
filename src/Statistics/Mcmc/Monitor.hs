{- |
Module      :  Statistics.Mcmc.Monitor
Description :  Monitor a Markov chain
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 21 14:35:11 2020.

-}

-- TODO: Probably name this Logger.

module Statistics.Mcmc.Monitor
  ( Out
  , Monitor
  , monitorFile
  , monitorScreen
  ) where

import qualified Data.ByteString as B
import Data.ByteString (ByteString)
import System.IO

-- | Either log to file or handle.
data Out = OutFile FilePath | OutHandle Handle

-- | Monitor a variable of the state space.
data Monitor a = Monitor
  {
    out       :: Out               -- ^ Log to file or handle ('stdout' can be used).
  , mShow     :: a -> ByteString   -- ^ Instruction about what to log.
  , frequency :: Int               -- ^ Logging frequency.
  }

-- | List of monitors.
data Supervision a = Supervision [Monitor a]

-- TODO: Write functions that act on 'Supervision', such as (open, print, close).

-- | Log to file.
monitorFile :: FilePath -> (a -> ByteString) -> Int -> Monitor a
monitorFile fp = Monitor (OutFile fp)

-- | Log to screen.
monitorScreen :: (a -> ByteString) -> Int -> Monitor a
monitorScreen = Monitor (OutHandle stdout)
