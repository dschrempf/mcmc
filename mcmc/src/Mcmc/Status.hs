-- Note: It is not necessary to add another type @b@ to store supplementary
-- information about the chain. The information can just be stored in @a@
-- equally well.

-- XXX: Status tuned exclusively to the Metropolis-Hastings algorithm. We should
-- abstract the algorithm from the chain. Maybe something like:
--
-- @
-- data Status a = Status { Chain a; Algorithm a}
-- @

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
  ( Cleaner (..),
    Status (..),
    status,
    cleanWith,
    saveWith,
    force,
    quiet,
    debug,
  )
where

import Data.Maybe
import Data.Time.Clock
import Mcmc.Item
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Trace
import Mcmc.Verbosity (Verbosity (..))
import Numeric.Log
import System.IO
import System.Random.MWC hiding (save)
import Prelude hiding (cycle)

-- | Clean the state periodically.
--
-- The prior and the likelihood will be updated after the cleaning process.
--
-- For long chains, successive numerical errors can accumulate such that the
-- state diverges from honoring specific required constraints. In these cases, a
-- 'Cleaner' can be used to ensure that the required constraints of the state
-- are honored. For example, the branches of an ultrametric phylogeny may
-- diverge slightly after successful many proposals such that the phylogeny is
-- not anymore ultrametric.
--
-- Please be aware that the Markov chain will not converge to the true posterior
-- distribution if the state is changed substantially! Only apply subtle changes
-- that are absolutely necessary to preserve the required properties of the
-- state such as specific numerical constraints.
data Cleaner a = Cleaner
  { -- | Clean every given number of iterations.
    clEvery :: Int,
    -- | Cleaning function. Executed before monitoring the state.
    clFunction :: a -> a
  }

-- | The 'Status' contains all information to run an MCMC chain. It is
-- constructed using the function 'status'.
--
-- The polymorphic type @a@ stores the state of the chain. It can also be used
-- to store auxiliary information.
data Status a = Status
  { -- MCMC related variables; saved.

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
    -- | For each 'Proposal', store the list of accepted (True) and rejected (False)
    -- proposals; for reasons of efficiency, the list is also stored in reverse
    -- order.
    acceptance :: Acceptance (Proposal a),
    -- | Number of burn in iterations; deactivate burn in with 'Nothing'.
    burnInIterations :: Maybe Int,
    -- | Auto tuning period (only during burn in); deactivate auto tuning with
    -- 'Nothing'.
    autoTuningPeriod :: Maybe Int,
    -- | Number of normal iterations excluding burn in. Note that auto tuning
    -- only happens during burn in.
    iterations :: Int,
    --
    -- Auxiliary variables; saved.

    -- | Overwrite output files? Default is 'False', change with 'force'.
    forceOverwrite :: Bool,
    -- | Save the chain with trace of given length at the end of the run?
    -- Default is no save ('Nothing'). Change with 'saveWith'.
    save :: Maybe Int,
    -- | Verbosity.
    verbosity :: Verbosity,
    -- | The random number generator.
    generator :: GenIO,
    --
    -- Auxiliary variables; not saved.

    -- | Starting time and starting iteration of chain; used to calculate
    -- run time and ETA.
    start :: Maybe (Int, UTCTime),
    -- | Handle to log file.
    logHandle :: Maybe Handle,
    --
    -- Auxiliary functions; not saved.

    -- | The prior function. The un-normalized posterior is the product of the
    -- prior and the likelihood.
    priorF :: a -> Log Double,
    -- | The likelihood function. The un-normalized posterior is the product of
    -- the prior and the likelihood.
    likelihoodF :: a -> Log Double,
    -- | Clean the state periodically.
    cleaner :: Maybe (Cleaner a),
    --
    -- Variables related to the algorithm; not saved.

    -- | A set of 'Proposal's form a 'Cycle'.
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
  -- | A list of 'Proposal's executed in forward order. The chain will be logged
  -- after each cycle.
  Cycle a ->
  -- | A 'Monitor' observing the chain.
  Monitor a ->
  -- | The initial state in the state space @a@.
  a ->
  -- | Number of burn in iterations; deactivate burn in with 'Nothing'.
  Maybe Int ->
  -- | Auto tuning period (only during burn in); deactivate auto tuning with
  -- 'Nothing'.
  Maybe Int ->
  -- | Number of normal iterations excluding burn in. Note that auto tuning only
  -- happens during burn in.
  Int ->
  -- | A source of randomness. For reproducible runs, make sure to use
  -- generators with the same, fixed seed.
  GenIO ->
  Status a
status n p l c m x mB mT nI g
  | isJust mT && isNothing mB = error "status: Auto tuning period given, but no burn in."
  | otherwise =
    Status
      n
      i
      0
      (singletonT i)
      (emptyA $ ccProposals c)
      mB
      mT
      nI
      False
      Nothing
      Info
      g
      Nothing
      Nothing
      p
      l
      Nothing
      c
      m
  where
    i = Item x (p x) (l x)

-- | Clean the state every given number of generations using the given function.
-- See 'Cleaner'.
cleanWith :: Cleaner a -> Status a -> Status a
cleanWith c s = s {cleaner = Just c}

-- | Save the Markov chain with trace of given length.
saveWith :: Int -> Status a -> Status a
saveWith n s = s {save = Just n}

-- | Overwrite existing files; it is not necessary to use 'force', when a chain
-- is continued.
force :: Status a -> Status a
force s = s {forceOverwrite = True}

-- | Do not print anything to standard output. Do not create log file. File
-- monitors and batch monitors are executed normally.
quiet :: Status a -> Status a
quiet s = s {verbosity = Quiet}

-- | Be verbose.
debug :: Status a -> Status a
debug s = s {verbosity = Debug}
