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

TODO: Maybe rename Chain to Trace; the Chain seems to be the list of States, and
the Trace includes the log-likelihoods and so on.

-}

module Statistics.Mcmc.Types
  ( -- * Static types
    --
    -- Types related to the status of an MCMC.
    State (..)
  , Mcmc (..)
  , Link (..)
  , Chain (..)
    -- * Dynamic types
    --
    -- Types related to moving between states.
  , Move (..)
  , Cycle (..)
  ) where

import Control.Monad.Primitive
import Numeric.Log
import System.Random.MWC

-- | A 'State' in the state space @a@.
newtype State a = State a
  deriving (Eq, Ord, Show, Read, Functor)

-- | An MCMC is defined by the current 'State' in the state space @a@, a
-- log-likelihood function, and a set of 'Move's forming a 'Cycle'. The wrapper
-- @m@ provides a source of randomness.
--
-- We could think of adding a set of moves, or something equivalent, but we
-- don't do that at the moment :).
data Mcmc a = Mcmc
  {
    -- ^ The current state.
    mcmcState :: State a
    -- ^ The log-likelihood function.
  , mcmcLlhf  :: State a -> Log Double
  }

-- | A 'Link' of the 'Chain'. For simplicity, we associate each passed 'State'
-- with the corresponding log-likelihood. A list of 'Link's is what makes a
-- 'Chain'.
data Link a = Link
  {
    -- ^ The current state.
    lnState   :: State a
    -- ^ The current log-likelihood.
  , lnLlh     :: Log Double
  }
  deriving (Eq, Ord, Show, Read)

-- | A 'Chain' passes through a list of 'State's with associated log-likelihoods
-- called 'Link's.
newtype Chain a = Chain [Link a]
  deriving (Show, Read)

instance Semigroup (Chain a) where
  (Chain l) <> (Chain r) = Chain (l <> r)

instance Monoid (Chain a) where
  mempty = Chain []

-- XXX: It may be advantageous to include these types into the Metropolis module.

-- | A 'Move' is an instruction about how the MCMC will traverse the 'State'
-- space. Essentially, it is a probability distribution conditioned on the
-- current 'State'.
--
-- We need to be able to sample a new state.
--
-- We need to know the probability of jumping forth, but also the probability of
-- jumping back. They are needed to calculate the Metropolis-Hastings ratio.
data Move m a = Move
  {
    -- ^ Instruction about sampling a new state.
    mvSample :: PrimMonad m => Gen (PrimState m) -> State a -> m (State a)
    -- ^ The log-likelihood of going from one state to another.
  , mvLlh :: State a -> State a -> Log Double
  }

-- | A collection of 'Move's form a 'Cycle'. The 'State' of the 'Chain' will be
-- logged after each 'Cycle'.
newtype Cycle m a = Cycle [Move m a]

instance Semigroup (Cycle m a) where
  (Cycle l) <> (Cycle r) = Cycle (l <> r)

instance Monoid (Cycle m a) where
  mempty = Cycle []

