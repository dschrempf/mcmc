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

-- TODO: Rename this module, Status is really bad ;).

-- TODO: Think about how to save and restore an MCMC run. It is easy to save and
-- restore the current state and likelihood (or the trace), but it seems
-- impossible to store all the moves and so on. This means, that one should
-- allow restart of a chain only with the same executable (which contains the
-- moves etc). See https://hackage.haskell.org/package/executable-hash.

-- TODO: Output ETA; time.

-- TODO: Output prior, likelihood, posterior.

module Statistics.Mcmc.Status
  ( Status (..)
  , getState
  , mcmc
  , reset
  , Mcmc
  , mcmcAutotune
  , mcmcOpenMonitors
  , mcmcExecMonitors
  , mcmcCloseMonitors
  ) where

import Prelude hiding (cycle)

import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict hiding (state)
import qualified Data.Map.Strict as M
import Data.Maybe
import Numeric.Log
import System.Random.MWC

import Statistics.Mcmc.Item
import Statistics.Mcmc.Monitor
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
    -- | A 'Monitor' observing the chain.
  , monitor       :: Monitor a
    -- | Number of completed cycles.
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

-- | Get current state of Markov chain.
getState :: Status a -> a
getState = state . item

-- | Initialize a Markov chain Monte Carlo run.
--
-- The 'Status' of a Markov chain includes information about the 'Move's, the
-- 'Trace', 'Acceptance' ratios, and more.
mcmc
  :: a                 -- ^ The initial state in the state space @a@.
  -> (a -> Log Double) -- ^ The un-normalized log-posterior function.
  -> Cycle a           -- ^ A list of 'Move's executed in forward order. The
                       -- chain will be logged after each cycle.
  -> Monitor a         -- ^ A 'Monitor' observing the chain.
  -> GenIO             -- ^ A source of randomness. For reproducible runs, make
                       -- sure to use a generator with the same seed.
  -> Status a          -- ^ The current 'Status' of the Markov chain.
mcmc x f c m = Status i f c m 0 (Trace [i]) (empty $ fromCycle c)
  where i   = Item x (f x)

-- | Reset a chain. Delete trace, acceptance ratios, and set the iteration to 0.
-- Used, for example, to reset a chain after burn in.
reset :: Status a -> Status a
reset s = s { iteration = 0, trace = Trace [i], acceptance = resetA a}
  where i = item s
        a = acceptance s

-- Tune the 'Move's in the 'Cycle' of the Markov chain 'Status'; check
-- acceptance ratio of the last n moves. Tuning has no effect on 'Move's that
-- cannot be tuned. See 'autotune'.
autotuneS :: Int -> Status a -> Status a
autotuneS n s = s {cycle = mapCycle tuneF (cycle s)}
  where
    ars  = acceptanceRatios n $ acceptance s
    tuneF m = fromMaybe m (autotune (ars M.! m) m)

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO

-- | Auto tune the 'Move's in the 'Cycle' of the chain. See 'autotune'.
mcmcAutotune :: Int -> Mcmc a ()
mcmcAutotune t = do
  liftIO $ putStrLn "-- Auto tune moves."
  modify' (autotuneS t)
  a <- gets acceptance
  c <- gets cycle
  liftIO $ putStr $ summarizeCycleA t a c

-- | Open the 'Monitor's of the chain. See 'msOpen'.
mcmcOpenMonitors :: Mcmc a ()
mcmcOpenMonitors = do
  liftIO $ putStrLn "-- Open monitors."
  s  <- get
  m  <- gets monitor
  m' <- liftIO $ mOpen m
  liftIO $ mHeader m'
  put s { monitor = m' }

-- | Execute the 'Monitor's of the chain. See 'msExec'.
mcmcExecMonitors :: Mcmc a ()
mcmcExecMonitors = do
  s <- get
  let i = iteration s
      x = getState s
      m = monitor s
  liftIO $ mExec i x m

-- | Close the 'Monitor's of the chain. See 'msClose'.
mcmcCloseMonitors :: Mcmc a ()
mcmcCloseMonitors = do
  liftIO $ putStrLn "-- Close monitors."
  m <- gets monitor
  liftIO $ mClose m
