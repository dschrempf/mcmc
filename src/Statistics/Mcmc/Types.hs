{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Statistics.Mcmc.Types
Description :  What is an MCMC?
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 18:01:15 2020.

TODO: Think about how to save and restore an MCMC run. It is easy to save and
restore the current state and likelihood (or the trace), but it seems impossible
to store all the moves and so on.

-}

module Statistics.Mcmc.Types
  (
    Item (..)
  , Trace (..)
  , prepend
  , Move (..)
  , Cycle (..)
  , Status (..)
  , start
  , Mcmc
  ) where

import Control.Monad.Trans.State.Strict
import Numeric.Log
import System.Random.MWC

-- | An 'State' of the 'Trace'. For simplicity, we associate each state 'S' with the
-- corresponding log-likelihood. A list of 'Item's is just a 'Trace'.
data Item a = Item
  {
    -- ^ The current state.
    -- lnState      :: S a
    state        :: a
    -- ^ The current log-posterior.
  , logPosterior :: Log Double
  }
  deriving (Eq, Ord, Show, Read)

-- | A 'Trace' passes through a list of states 'S' with associated log-likelihoods
-- called 'Item's.
newtype Trace a = Trace [Item a]
  deriving (Show, Read)

instance Semigroup (Trace a) where
  (Trace l) <> (Trace r) = Trace (l <> r)

instance Monoid (Trace a) where
  mempty = Trace []

-- Prepend an 'Item' to a 'Trace'.
{-# INLINE prepend #-}
prepend :: Item a -> Trace a -> Trace a
prepend x (Trace xs) = Trace (x:xs)

-- | A 'Move' is an instruction about how the MCMC will traverse the state
-- space. Essentially, it is a probability distribution conditioned on the
-- current state 'S'.
--
-- We need to be able to sample a new state.
--
-- We need to know the probability of jumping forth, but also the probability of
-- jumping back. They are needed to calculate the Metropolis-Hastings ratio.
data Move a = Move
  {
    -- TODO: Don't use IO explicitly here.
    --
    -- ^ Instruction about sampling a new state.
    mvSample   :: a -> GenIO -> IO a
    -- ^ The log-probability of going from one state to another.
  , mvLogProb  :: a -> a -> Log Double
  }

-- | A collection of 'Move's form a 'Cycle'. The state 'S' of the 'Trace' will
-- be logged after each 'Cycle'.
newtype Cycle a = Cycle [Move a]

instance Semigroup (Cycle a) where
  (Cycle l) <> (Cycle r) = Cycle (l <> r)

instance Monoid (Cycle a) where
  mempty = Cycle []

-- | The 'Status' of an MCMC run.
data Status a = Status
  {
    -- ^ The current or initial 'Item' of the chain combines the current state
    -- 'S' and the current log-likelihood.
    item          :: Item a
    -- ^ The unnormalized log-posterior log-likelihood function.
  , logPosteriorF :: a -> Log Double
    -- ^ A set of 'Move's form a 'Cycle'.
  , moves         :: Cycle a
    -- ^ The 'Trace' of the Markov chain in reverse order, newest first.
  , trace         :: Trace a
    -- ^ The random number generator.
  , generator     :: GenIO
  }

-- | Initialize a Markov chain Monte Carlo run with a specified seed.
start :: a -> (a -> Log Double) -> Cycle a -> GenIO -> Status a
start x f c = Status i f c (Trace [i])
  where i = Item x (f x)

-- | An Mcmc state transformer.
--
-- TODO: Improve documentation.
type Mcmc a = StateT (Status a) IO
