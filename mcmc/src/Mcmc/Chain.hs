-- |
-- Module      :  Mcmc.Chain
-- Description :  What is an MCMC?
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 18:01:15 2020.
module Mcmc.Chain
  ( Chain (..),
    chain,
    noData,
    Cleaner (..),
    cleanWith,
  )
where

-- Note: It is not necessary to add another type @b@ to store supplementary
-- information about the chain. The information can just be stored in @a@
-- equally well.

-- XXX: Status tuned exclusively to the Metropolis-Hastings algorithm. We should
-- abstract the algorithm from the chain. Maybe something like:
--
-- @
-- data Status a = Status { Chain a; Algorithm a}
-- @

-- TODO: REFACTOR. Check documentation.

import Data.Maybe
import Data.Time.Clock
import Mcmc.Item
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Trace
import Numeric.Log
import System.Random.MWC hiding (save)
import Prelude hiding (cycle)

-- | The 'Chain' contains all information to run a Markov chain Monte Carlo
-- sampler. A 'Chain' is constructed using the function 'chain'.
--
-- The state of the chain has type @a@. If necessary, the type @a@ can also be
-- used to store auxiliary information.
--
-- The 'Chain' stores information about current state ('Mcmc.Item.Item') and
-- iteration, the history of the chain ('Mcmc.Trace.Trace'), the
-- 'Acceptance' ratios, and the random number generator.
--
-- Further, the 'Chain' includes auxiliary variables and functions such as
-- the prior and likelihood functions, instructions to move around the state
-- space (see above) and to monitor the MCMC run.
--
-- The 'Mcmc.Environment.Environment' of the chain is excluded.
data Chain a = Chain
  { -- Variables; saved.

    -- | The name of the chain; used as file prefix.
    chainName :: String,
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
    -- | The random number generator.
    generator :: GenIO,
    --
    -- Variables and functions; not saved.

    -- | Starting time and starting iteration of chain; used to calculate
    -- run time and ETA.
    start :: Maybe (Int, UTCTime),
    -- | The prior function. The un-normalized posterior is the product of the
    -- prior and the likelihood.
    priorF :: a -> Log Double,
    -- | The likelihood function. The un-normalized posterior is the product of
    -- the prior and the likelihood.
    likelihoodF :: a -> Log Double,
    -- | Clean the state periodically.
    cleaner :: Maybe (Cleaner a),
    -- | A set of 'Proposal's form a 'Cycle'.
    --
    -- TODO: Should we move the cycle to a dedicated @Algorithm@ type?
    cycle :: Cycle a,
    -- | A 'Monitor' observing the chain.
    monitor :: Monitor a
  }

-- | Initialize a Markov chain.
chain ::
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
  -- | A source of randomness. For reproducible runs, make sure to use
  -- generators with the same, fixed seed.
  GenIO ->
  Chain a
chain n p l c m x g =
  Chain
    n
    i
    0
    (singletonT i)
    (emptyA $ ccProposals c)
    g
    Nothing
    p
    l
    Nothing
    c
    m
  where
    i = Item x (p x) (l x)

-- | Set the likelihood function to 1.0. Useful for debugging and testing.
noData :: Chain a -> Chain a
noData x = x {likelihoodF = const 1.0}

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

-- | Clean the state every given number of generations using the given function.
-- See 'Cleaner'.
cleanWith :: Cleaner a -> Chain a -> Chain a
cleanWith c x = x {cleaner = Just c}
