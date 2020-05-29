{- |
Module      :  Statistics.Mcmc.Status
Description :  What is an MCMC?
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 18:01:15 2020.

-}

-- TODO: Think about how to save and restore an MCMC run. It is easy to save and
-- restore the current state and likelihood (or the trace), but it seems
-- impossible to store all the moves and so on. This means, that one should
-- allow restart of a chain only with the same executable (which contains the
-- moves etc). See https://hackage.haskell.org/package/executable-hash.

module Statistics.Mcmc.Status
  ( Status (..)
  , status
  ) where

import Prelude hiding (cycle)

import Numeric.Log
import System.Random.MWC

import Statistics.Mcmc.Item
import Statistics.Mcmc.Monitor
import Statistics.Mcmc.Move
import Statistics.Mcmc.Trace

-- TODO: Add possibility to store supplementary information about the chain.

-- | The 'Status' of an MCMC run; see 'status' for creation.
data Status a = Status
  {
    -- | The current 'Item' of the chain combines the current state and the
    -- current log-likelihood.
    item          :: Item a
    -- | The un-normalized log-posterior function. The log-posterior is the sum
    -- of the log-prior and the log-likelihood.
  , logPosteriorF :: a -> Log Double
    -- | A set of 'Move's form a 'Cycle'.
  , cycle         :: Cycle a
    -- | A 'Monitor' observing the chain.
  , monitor       :: Monitor a
    -- | The iteration is the number of completed cycles.
  , iteration     :: Int
    -- | The 'Trace' of the Markov chain in reverse order, the most recent
    -- 'Item' is at the head of the list.
  , trace         :: Trace a
    -- | For each 'Move', store the list of accepted (True) and rejected (False)
    -- proposals; for reasons of efficiency, the list is also stored in reverse
    -- order.
  , acceptance    :: Acceptance (Move a)
    -- | The random number generator.
  , generator     :: GenIO
  }

-- | Initialize the status of a Markov chain Monte Carlo run.
--
-- The 'Status' of a Markov chain includes information about the 'Move's, the
-- 'Trace', 'Acceptance' ratios, and more.
status
  :: a                 -- ^ The initial state in the state space @a@.
  -> (a -> Log Double) -- ^ The un-normalized log-posterior function.
  -> Cycle a           -- ^ A list of 'Move's executed in forward order. The
                       -- chain will be logged after each cycle.
  -> Monitor a         -- ^ A 'Monitor' observing the chain.
  -> GenIO             -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> Status a          -- ^ The current 'Status' of the Markov chain.
status x f c m = Status i f c m 0 (Trace [i]) (empty $ fromCycle c)
  where i   = Item x (f x)

-- -- | Get current state of Markov chain.
-- getState :: Status a -> a
-- getState = state . item

-- -- | Reset a chain. Delete trace, acceptance ratios, and set the iteration to 0.
-- -- Used, for example, to reset a chain after burn in.
-- reset :: Status a -> Status a
-- reset s = s { iteration = 0, trace = Trace [i], acceptance = resetA a}
--   where i = item s
--         a = acceptance s
