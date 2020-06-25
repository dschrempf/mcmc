-- XXX: Add possibility to store supplementary information about the chain.
--
-- Maybe something like Trace b; and give a function a -> b to extract
-- supplementary info.

-- XXX: Status tuned exclusively to the Metropolis-Hastings algorithm. We
-- should abstract the algorithm from the chain. For example,
--
-- @
-- data Status a b = Status { Chain a; Algorithm a b}
-- @
--
-- where a described the state space and b the auxiliary information of the
-- algorithm. This would also solve the above problem, for example in terms of
-- the Hamiltonian algorithm

-- |
-- Module      :  Mcmc.Status
-- Description :  What is an MCMC?
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 18:01:15 2020.
module Mcmc.Status
  ( Status (..),
    status,
    noSave,
  )
where

import Data.Time.Clock
import Mcmc.Item
import Mcmc.Monitor
import Mcmc.Move
import Mcmc.Trace
import Numeric.Log
import System.Random.MWC hiding (save)
import Prelude hiding (cycle)

-- | The 'Status' contains all information to run an MCMC chain. It is
-- constructed using the function 'status'.
data Status a
  = Status
      { -- Variables saved to disc.

        -- | The name of the MCMC chain; used as file prefix.
        name :: String,
        -- | The current 'Item' of the chain combines the current state and the
        -- current likelihood.
        item :: Item a,
        -- | The iteration is the number of completed cycles.
        iteration :: Int,
        -- | The 'Trace' of the Markov chain in reverse order, the most recent
        -- 'Item' is at the head of the list.
        trace :: Trace a,
        -- | For each 'Move', store the list of accepted (True) and rejected (False)
        -- proposals; for reasons of efficiency, the list is also stored in reverse
        -- order.
        acceptance :: Acceptance (Move a),
        -- | Number of burn in iterations; deactivate burn in with 'Nothing'.
        burnInIterations :: Maybe Int,
        -- | Auto tuning period (only during burn in); deactivate auto tuning with
        -- 'Nothing'.
        autoTuningPeriod :: Maybe Int,
        -- | Number of normal iterations excluding burn in. Note that auto tuning
        -- only happens during burn in.
        iterations :: Int,
        -- | Starting time and starting iteration of chain; used to calculate
        -- run time and ETA.
        start :: Maybe (Int, UTCTime),
        -- | Save the chain at the end of the run? Defaults to 'True'.
        save :: Bool,
        -- | The random number generator.
        generator :: GenIO,
        -- Auxiliary functions.

        -- | The prior function. The un-normalized posterior is the product of the
        -- prior and the likelihood.
        priorF :: a -> Log Double,
        -- | The likelihood function. The un-normalized posterior is the product of
        -- the prior and the likelihood.
        likelihoodF :: a -> Log Double,
        -- Variables related to the algorithm.

        -- | A set of 'Move's form a 'Cycle'.
        cycle :: Cycle a,
        -- | A 'Monitor' observing the chain.
        monitor :: Monitor a
      }

-- | Initialize the 'Status' of a Markov chain Monte Carlo run.
status ::
  -- | Name of the Markov chain; used as file prefix.
  String ->
  -- | The prior function.
  (a -> Log Double) ->
  -- | The likelihood function.
  (a -> Log Double) ->
  -- | A list of 'Move's executed in forward order. The
  -- chain will be logged after each cycle.
  Cycle a ->
  -- | A 'Monitor' observing the chain.
  Monitor a ->
  -- | The initial state in the state space @a@.
  a ->
  -- | Number of burn in iterations; deactivate burn in with 'Nothing'.
  Maybe Int ->
  -- | Auto tuning period (only during burn in); deactivate
  -- auto tuning with 'Nothing'.
  Maybe Int ->
  -- | Number of normal iterations excluding burn in. Note
  -- that auto tuning only happens during burn in.
  Int ->
  -- | A source of randomness. For reproducible runs, make
  -- sure to use a generator with the same seed.
  GenIO ->
  Status a
status n p l c m x mB mT nI g =
  Status
    n
    i
    0
    (singletonT i)
    (emptyA $ ccMoves c)
    mB
    mT
    nI
    Nothing
    True
    g
    p
    l
    c
    m
  where
    i = Item x (p x) (l x)

-- | Do not save the Markov chain at the end.
noSave :: Status a -> Status a
noSave s = s {save = False}
