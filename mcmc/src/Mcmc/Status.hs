{- |
Module      :  Mcmc.Status
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
--
-- Idea: Separate Status into information stored on disk, and retrieved upon
-- continuing an analysis (item, trace, acceptance, generator, iteration;
-- something else?). Possibly limit the length of the trace to the maximum batch
-- size.
--
-- The moves, the posterior, and so on, have to be provided again for the next
-- analysis. Recompute and check the posterior for the last state because the
-- posterior function may have changed. Of course, we cannot test for the same
-- function, but having the same posterior at the last state is already a good
-- indicator.
--
-- Results: I tried this, and it is slow, too slow. The reason is that the
-- compound state, which I called McmcData, has to be unpacked and repacked all
-- the time.
--
-- The next idea is: write a function of type 'Status a -> Save a'; and a
-- continue function, something like 'Save a -> auxiliary information -> Status
-- a', where 'Save a' is a data type combining all the information to be stored.

module Mcmc.Status
  ( Status(..)
  , status
  )
where

import           Prelude                 hiding ( cycle )

import           Data.Time.Clock
import           Numeric.Log
import           System.Random.MWC

import           Mcmc.Item
import           Mcmc.Monitor
import           Mcmc.Move
import           Mcmc.Trace

-- TODO: Add possibility to store supplementary information about the chain.

-- | The 'Status' of an MCMC run; see 'status' for creation.
data Status a = Status
  {
    -- Information stored at each generation.
    -- | The current 'Item' of the chain combines the current state and the
    -- current log-likelihood.
    item            :: Item a
    -- | The iteration is the number of completed cycles.
  , iteration       :: Int
    -- | For each 'Move', store the list of accepted (True) and rejected (False)
    -- proposals; for reasons of efficiency, the list is also stored in reverse
    -- order.
  , acceptance      :: Acceptance (Move a)
    -- | The 'Trace' of the Markov chain in reverse order, the most recent
    -- 'Item' is at the head of the list.
  , trace           :: Trace a
    -- | The random number generator.
  , generator       :: GenIO

    -- Auxiliary information about the chain, which has to be provided again
    -- upon restart or continuation.
    -- | The log-prior function. The un-normalized log-posterior is the sum of
    -- the log-prior and the log-likelihood.
  , logPriorF       :: a -> Log Double
    -- | The log-likelihood function. The un-normalized log-posterior is the sum
    -- of the log-prior and the log-likelihood.
  , logLikelihoodF  :: a -> Log Double
    -- | A set of 'Move's form a 'Cycle'.
  , cycle           :: Cycle a
    -- | A 'Monitor' observing the chain.
  , monitor         :: Monitor a
    -- | Total number of iterations.
  , totalIterations :: Maybe Int
    -- | Starting time of chain; used to calculate run time and ETA.
  , starttime       :: Maybe UTCTime
  }

-- | Initialize the status of a Markov chain Monte Carlo run.
--
-- The 'Status' of a Markov chain includes information about the 'Move's, the
-- 'Trace', 'Acceptance' ratios, and more.
status
  :: a                 -- ^ The initial state in the state space @a@.
  -> (a -> Log Double) -- ^ The log-prior function.
  -> (a -> Log Double) -- ^ The log-likelihood function.
  -> Cycle a           -- ^ A list of 'Move's executed in forward order. The
                       -- chain will be logged after each cycle.
  -> Monitor a         -- ^ A 'Monitor' observing the chain.
  -> GenIO             -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> Status a          -- ^ The current 'Status' of the Markov chain.
status x p l c m g = Status i
                            0
                            (empty $ fromCycle c)
                            (Trace [i])
                            g
                            p
                            l
                            c
                            m
                            Nothing
                            Nothing
  where i = Item x (p x) (l x)
