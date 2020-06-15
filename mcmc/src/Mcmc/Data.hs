{- |
Module      :  Mcmc.Data
Description :  What is an MCMC?
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 18:01:15 2020.

-}

-- TODO: Continue an MCMC run. It is easy to save and restore the current state
-- and likelihood (or the trace), but it is not feasible to store all the moves
-- and so on, so they have to be provided again ('Spec' data type).

-- TODO: upon continuation: recompute and check the posterior for the last state
-- because the posterior function may have changed. Of course, we cannot test
-- for the same function, but having the same posterior at the last state is
-- already a good indicator.

-- TODO: Possibly limit the length of the trace to the maximum batch size.

module Mcmc.Data
  ( McmcData(..)
  , mcmc
  , Spec(..)
  , spec
  , Status(..)
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

-- | All we need to run a chain combined in one data type. 'Status' and 'Spec'
-- are separated such that storable information can be written to disc and
-- restored to continue a Markov chain run.
data McmcData a = McmcData
  {
    mcmcStatus :: Status a
  , mcmcSpec   :: Spec a
  }

-- | All we need to run a chain combined in one data type.
mcmc
  :: (a -> Log Double) -- ^ The log-prior function.
  -> (a -> Log Double) -- ^ The log-likelihood function.
  -> Cycle a           -- ^ A list of 'Move's executed in forward order. The
                       -- chain will be logged after each cycle.
  -> Monitor a         -- ^ A 'Monitor' observing the chain.
  -> a                 -- ^ The initial state in the state space @a@.
  -> GenIO             -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> McmcData a
mcmc p l c m x g = McmcData st sp
  where sp = spec p l c m
        st = status sp x g

-- | The 'Spec'ification of an MCMC run. The 'Spec' includes the log prior and
-- log likelihood functions, instructions to move around the state space and to
-- monitor the MCMC run, as well as some auxiliary information. See 'spec' for
-- creation.
data Spec a = Spec
  {
    -- | The log-prior function. The un-normalized log-posterior is the sum of
    -- the log-prior and the log-likelihood.
    logPriorF       :: a -> Log Double
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

-- | Initialize the 'Spec'ification of a Markov chain Monte Carlo run.
--
-- The 'Spec' of a Markov chain includes instructions about how to calculate the
-- prior, and the likelihood, about how to 'Move' around the state space, and
-- what should be monitored, as well as other auxiliary information.
spec
  :: (a -> Log Double) -- ^ The log-prior function.
  -> (a -> Log Double) -- ^ The log-likelihood function.
  -> Cycle a           -- ^ A list of 'Move's executed in forward order. The
                       -- chain will be logged after each cycle.
  -> Monitor a         -- ^ A 'Monitor' observing the chain.
  -> Spec a
spec p l c m = Spec p l c m Nothing Nothing

-- TODO: Add possibility to store supplementary information about the chain.
-- Should I put the information into Spec or Status?
--
-- Maybe something like Trace b; and give a function a -> b to extract
-- supplementary info.

-- | The 'Status' of an MCMC run. Information that can be saved to disc, so that
-- a chain can be continued. See 'status' for creation.
data Status a = Status
  {
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
  }

-- TODO.
--
-- instance (JSON a) => JSON (Status a)

-- | Initialize the 'Status' of a Markov chain Monte Carlo run.
--
-- The 'Status' of a Markov chain includes information about current state
-- ('Item'), the 'Trace', the 'Acceptance' ratios, and the random number
-- generator.
status
  :: Spec a
  -> a                 -- ^ The initial state in the state space @a@.
  -> GenIO             -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> Status a          -- ^ The current 'Status' of the Markov chain.
status sp x = Status i 0 (Trace [i]) (empty $ fromCycle c)
  where
    p = logPriorF sp
    l = logLikelihoodF sp
    c = cycle sp
    i = Item x (p x) (l x)
