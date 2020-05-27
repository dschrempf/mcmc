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
-- restore the current state and likelihood (or the trace), but it seems impossible
-- to store all the moves and so on.

-- TODO: Define monitors like in RevBayes.

-- TODO: Proper output. RevBayes does the following:
-- @
--   Running burn-in phase of Monte Carlo sampler for 10000 iterations.
--   This simulation runs 1 independent replicate.
--   The MCMCMC simulator runs 1 cold chain and 3 heated chains.
--   The simulator uses 42 different moves in a random move schedule with 163 moves per iteration
-- @

module Statistics.Mcmc.Status
  ( Status (..)
  , mcmc
  , reset
  , autotuneS
  , Mcmc
  ) where

import Prelude hiding (cycle)

import Control.Monad.Trans.State.Strict
import qualified Data.Map.Strict as M
import Data.Maybe
import Numeric.Log
import System.Random.MWC

import Statistics.Mcmc.Acceptance
import Statistics.Mcmc.Item
import Statistics.Mcmc.Move
import Statistics.Mcmc.Trace

-- TODO: Add possibility to store supplementary information about the chain.

-- | The 'Status' of an MCMC run.
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
    -- | Number of completed cycles.
  , iteration     :: Int
    -- | The 'Trace' of the Markov chain in reverse order, the most recent state
    -- is at the head of the list.
  , trace         :: Trace a
    -- | For each 'Move', store the list of accepted (True) and rejected (False)
    -- proposals; for reasons of efficiency, the list is also stored in reverse
    -- order.
  , acceptance    :: Acceptance (Move a)
    -- | The random number generator.
  , generator     :: GenIO
  }

-- | Initialize a Markov chain Monte Carlo run.
--
-- The 'Status' of a Markov chain includes information about the 'Move's, the
-- 'Trace', 'Acceptance' ratios, and more.
mcmc
  :: a -- ^ The initial state in the state space @a@.
  -> (a -> Log Double) -- ^ The un-normalized log-posterior function.
  -> Cycle a -- ^ A list of 'Move's executed in forward order. The chain will
               -- be logged after each cycle.
  -> GenIO -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> Status a -- ^ The current 'Status' of the Markov chain.
mcmc x f c = Status i f c 0 (Trace [i]) (empty $ moves c)
  where i   = Item x (f x)

-- | Reset a chain. Delete trace, acceptance ratios, and set the iteration to 0.
-- Used, for example, to reset a chain after burn in.
reset :: Status a -> Status a
reset s = s { iteration = 0, trace = Trace [i], acceptance = resetA a}
  where i = item s
        a = acceptance s

-- | Tune the 'Move's in the 'Cycle' of the Markov chain 'Status'; check
-- acceptance ratio of the last n moves. Tuning has no effect on 'Move's that
-- cannot be tuned. See 'autotune'.
autotuneS :: Int -> Status a -> Status a
autotuneS n s = s {cycle = mapC tuneF (cycle s)}
  where
    ars  = acceptanceRatios n $ acceptance s
    tuneF m = fromMaybe m (autotune (ars M.! m) m)

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO

