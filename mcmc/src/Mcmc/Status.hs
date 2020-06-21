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

-- TODO: Add possibility to store supplementary information about the chain.
--
-- Maybe something like Trace b; and give a function a -> b to extract
-- supplementary info.

-- TODO: Status tuned exclusively to the Metropolis-Hastings algorithm. We
-- should abstract the algorithm from the chain. For example,
--
-- @
-- data Status a b = Status { Chain a; Algorithm a b}
-- @
--
-- where a described the state space and b the auxiliary information of the
-- algorithm. This would also solve the above problem, for example in terms of
-- the Hamiltonian algorithm

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
    -- | The name of the MCMC chain; used as file prefix.
    name             :: String
    -- | The current 'Item' of the chain combines the current state and the
    -- current log-likelihood.
  , item             :: Item a
    -- | The iteration is the number of completed cycles.
  , iteration        :: Int
    -- | The 'Trace' of the Markov chain in reverse order, the most recent
    -- 'Item' is at the head of the list.
  , trace            :: Trace a
    -- | For each 'Move', store the list of accepted (True) and rejected (False)
    -- proposals; for reasons of efficiency, the list is also stored in reverse
    -- order.
  , acceptance       :: Acceptance (Move a)
    -- | Number of burn in iterations; deactivate burn in with 'Nothing'.
  , burnInIterations :: Maybe Int
    -- | Auto tuning period (only during burn in); deactivate auto tuning with
    -- 'Nothing'.
  , autoTuningPeriod :: Maybe Int
    -- | Number of normal iterations excluding burn in. Note that auto tuning
    -- only happens during burn in.
  , iterations       :: Int
    -- | Starting time of chain; used to calculate run time and ETA.
  , starttime        :: Maybe UTCTime
    -- | time of last save.
  , savetime         :: Maybe UTCTime
    -- | The random number generator.
  , generator        :: GenIO

    -- Auxiliary functions.
    -- | The log-prior function. The un-normalized log-posterior is the sum of
    -- the log-prior and the log-likelihood.
  , logPriorF        :: a -> Log Double
    -- | The log-likelihood function. The un-normalized log-posterior is the sum
    -- of the log-prior and the log-likelihood.
  , logLikelihoodF   :: a -> Log Double

    -- Variables related to the algorithm.
    -- | A set of 'Move's form a 'Cycle'.
  , cycle            :: Cycle a
    -- | A 'Monitor' observing the chain.
  , monitor          :: Monitor a
  }

-- | Initialize the 'Status' of a Markov chain Monte Carlo run.
status
  :: String            -- ^ Name of the Markov chain; used as file prefix.
  -> (a -> Log Double) -- ^ The log-prior function.
  -> (a -> Log Double) -- ^ The log-likelihood function.
  -> Cycle a           -- ^ A list of 'Move's executed in forward order. The
                       -- chain will be logged after each cycle.
  -> Monitor a         -- ^ A 'Monitor' observing the chain.
  -> a                 -- ^ The initial state in the state space @a@.
  -> Maybe Int         -- ^ Number of burn in iterations; deactivate burn in with 'Nothing'.
  -> Maybe Int         -- ^ Auto tuning period (only during burn in); deactivate
                       -- auto tuning with 'Nothing'.
  -> Int               -- ^ Number of normal iterations excluding burn in. Note
                       -- that auto tuning only happens during burn in.
  -> GenIO             -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> Status a
status n p l c m x mB mT nI g = Status n
                                       i
                                       0
                                       (singletonT i)
                                       (emptyA $ fromCycle c)
                                       mB
                                       mT
                                       nI
                                       Nothing
                                       Nothing
                                       g
                                       p
                                       l
                                       c
                                       m
  where i = Item x (p x) (l x)