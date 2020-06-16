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

-- TODO: Possibly limit the length of the trace to the maximum batch size.

-- TODO: Add possibility to store supplementary information about the chain.
--
-- Maybe something like Trace b; and give a function a -> b to extract
-- supplementary info.

module Mcmc.Status
  (
    Status(..)
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

-- | The 'Status' of an MCMC run. All we need to run a chain combined in one
-- data type. See 'status' for creation.
--
-- The 'Status' of a Markov chain includes information about current state
-- ('Item') and iteration, the history of the chain ('Trace'), the 'Acceptance'
-- ratios, and the random number generator.
--
-- Further, the 'Status' includes auxiliary variables and functions such as the
-- log prior and log likelihood functions, instructions to move around the state
-- space and to monitor the MCMC run, as well as some auxiliary information.
data Status a = Status
  {
    -- Variables saved to disc.
    -- | The current 'Item' of the chain combines the current state and the
    -- current log-likelihood.
    item            :: Item a
    -- | The iteration is the number of completed cycles.
  , iteration       :: Int
    -- | The 'Trace' of the Markov chain in reverse order, the most recent
    -- 'Item' is at the head of the list.
  , trace           :: Trace a
    -- | For each 'Move', store the list of accepted (True) and rejected (False)
    -- proposals; for reasons of efficiency, the list is also stored in reverse
    -- order.
  , acceptance      :: Acceptance (Move a)
    -- | The random number generator.
  , generator       :: GenIO

    -- Auxiliary variables and functions.
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
    -- | Starting time of chain; used to calculate run time and ETA.
  , starttime       :: Maybe UTCTime
    -- | Total number of iterations.
  , totalIterations :: Maybe Int
  }

-- | Initialize the 'Status' of a Markov chain Monte Carlo run.
status
  :: (a -> Log Double) -- ^ The log-prior function.
  -> (a -> Log Double) -- ^ The log-likelihood function.
  -> Cycle a           -- ^ A list of 'Move's executed in forward order. The
                       -- chain will be logged after each cycle.
  -> Monitor a         -- ^ A 'Monitor' observing the chain.
  -> a                 -- ^ The initial state in the state space @a@.
  -> GenIO             -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> Status a
status p l c m x g = Status i 0 (Trace [i]) (empty $ fromCycle c) g p l c m Nothing Nothing
  where
    i = Item x (p x) (l x)
