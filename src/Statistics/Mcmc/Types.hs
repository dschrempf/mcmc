{-# LANGUAGE DeriveFunctor #-}
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
  ( S (..)
  , Item (..)
  , Trace (..)
  , Move (..)
  , Cycle (..)
  , LogLikelihoodFunction
  , Status (..)
  -- , Mcmc (..)
  ) where

import Control.Monad.Primitive
import Numeric.Log
import System.Random.MWC

-- | A state 'S' in the state space @a@.
newtype S a = S a
  deriving (Eq, Ord, Show, Read, Functor)

-- | An 'Item' of the 'Trace'. For simplicity, we associate each state 'S' with the
-- corresponding log-likelihood. A list of 'Item's is just a 'Trace'.
data Item a = Item
  {
    -- ^ The current state.
    lnState   :: S a
    -- ^ The current log-likelihood.
  , lnLlh     :: Log Double
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

-- | A 'Move' is an instruction about how the MCMC will traverse the state
-- space. Essentially, it is a probability distribution conditioned on the
-- current state 'S'.
--
-- We need to be able to sample a new state.
--
-- We need to know the probability of jumping forth, but also the probability of
-- jumping back. They are needed to calculate the Metropolis-Hastings ratio.
data Move m a = Move
  {
    -- ^ Instruction about sampling a new state.
    mvSample :: PrimMonad m => S a -> Gen (PrimState m) -> m (S a)
    -- ^ The log-likelihood of going from one state to another.
  , mvLlh :: S a -> S a -> Log Double
  }

-- | A collection of 'Move's form a 'Cycle'. The state 'S' of the 'Trace' will
-- be logged after each 'Cycle'.
newtype Cycle m a = Cycle [Move m a]

instance Semigroup (Cycle m a) where
  (Cycle l) <> (Cycle r) = Cycle (l <> r)

instance Monoid (Cycle m a) where
  mempty = Cycle []

-- | The log-likelihood function maps a state 'S' to a log-likelihood.
type LogLikelihoodFunction a = S a -> Log Double

-- | The 'Status' of an MCMC run is defined by
--   - the current state in the state space @a@ together with the
--     log-likelihood,
--   - a log-likelihood function,
--   - a set of 'Move's forming a 'Cycle', and
--   - the 'Trace'.
--
-- The wrapper @m@ provides a source of randomness.
data Status m a = Status
  {
    -- ^ The current 'Item' of the chain combines the current state 'S' and
    -- log-likelihood.
    mcmcState :: Item a
    -- ^ The log-likelihood function.
  , mcmcLlhf  :: LogLikelihoodFunction a
    -- ^ A set of 'Move's form a 'Cycle'.
  , mcmcCycle :: Cycle m a
    -- ^ The 'Trace' of the Markov chain in reverse order, newest first.
  , mcmcTrace :: Trace a
  }

-- -- | An Mcmc state transformer.
-- --
-- -- TODO: Improve documentation.
-- newtype Mcmc m a = StateT (Status m a)
