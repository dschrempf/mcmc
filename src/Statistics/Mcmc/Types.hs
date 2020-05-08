{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE DeriveFunctor #-}

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
  , prependT
  , prependA
  , Move (..)
  , Cycle (..)
  , Status (..)
  , start
  , Mcmc
  ) where

import Control.Monad.Trans.State.Strict
import Numeric.Log
import System.Random.MWC

-- | An 'Item', or link of the Markov chain. For reasons of computational
-- efficiency, each state is associated with the corresponding log-likelihood. A
-- list of 'Item's is called the 'Trace' of the Markov chain.
data Item a = Item
  {
    -- | The current state in the state space @a@.
    state        :: a
    -- | The current log-posterior.
  , logPosterior :: Log Double
  }
  deriving (Eq, Ord, Show, Read)

-- | A 'Trace' passes through a list of states with associated log-likelihoods
-- which are called 'Item's. New 'Item's are prepended, and the path of the
-- Markov chain is stored in reversed order.
newtype Trace a = Trace {fromTrace :: [Item a] }
  deriving (Show, Read)

instance Semigroup (Trace a) where
  (Trace l) <> (Trace r) = Trace (l <> r)

instance Monoid (Trace a) where
  mempty = Trace []

-- | Prepend an 'Item' to a 'Trace'.
{-# INLINE prependT #-}
prependT :: Item a -> Trace a -> Trace a
prependT x (Trace xs) = Trace (x:xs)

-- TODO: This is highly confusing. Make this more clear. Probably use a data
-- type.

-- | Prepend a list of accepted / rejected tries to the list of
-- accepted / rejected tries.
{-# INLINE prependA #-}
prependA :: [Bool] -> [[Bool]] -> [[Bool]]
prependA as ass = [ x:xs | (x, xs) <- zip as ass]

-- | A 'Move' is an instruction about how the Markov chain will traverse the
-- state space @a@. Essentially, it is a probability density conditioned on the
-- current state.
--
-- We need to know the probability density of jumping forth, but also the
-- probability density of jumping back. They are needed to calculate the
-- Metropolis-Hastings ratio.
data Move a = Move
  {
    -- | The name of the move.
    mvName       :: String
    -- | Instruction about randomly moving to a new state.
  , mvSample     :: a -> GenIO -> IO a
    -- | The log-density of going from one state to another.
  , mvLogDensity :: a -> a -> Log Double
  }

-- | A collection of 'Move's form a 'Cycle'. The state of the 'Trace' will be
-- logged only after each 'Cycle'.
--
-- XXX: At the moment, 'Move's are executed in forward order as they appear in
-- the 'Cycle'. One could think of associating weighs with moves, and execute
-- moves according to their relative weights in the Cycle --- a common practice
-- in MCMC software packages.
newtype Cycle a = Cycle { fromCycle :: [Move a] }

instance Semigroup (Cycle a) where
  (Cycle l) <> (Cycle r) = Cycle (l <> r)

instance Monoid (Cycle a) where
  mempty = Cycle []

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
  , moves         :: Cycle a
    -- | The 'Trace' of the Markov chain in reverse order, the most recent state
    -- is at the head of the list.
  , trace         :: Trace a
    -- TODO: See 'prependA'.
    --
    -- | Log of accepted (True) and rejected (False) 'Move's per 'Cycle'; also
    -- stored in reverse order. @acceptance status !! i@ is the list of
    -- accepted/rejected tries of the @i@th 'Move' in the 'Cycle'.
  , acceptance    :: [[Bool]]
    -- | The random number generator.
  , generator     :: GenIO
  }

-- | Initialize a Markov chain Monte Carlo run.
start
  :: a -- ^ The initial state in the state space @a@.
  -> (a -> Log Double) -- ^ The un-normalized log-posterior function.
  -> Cycle a -- ^ A list of 'Move's executed in forward order. The chain will be
             -- logged after each cycle.
  -> GenIO -- ^ A source of randomness. For reproducible runs, make sure to use
           -- a generator with the same seed.
  -> Status a -- ^ Return the current 'Status' of the Markov chain. Use, for
              -- example, 'Statistics.Mcmc.Metropolis.mh' to run the Markov
              -- chain.
start x f c = Status i f c (Trace [i]) (replicate nMoves [])
  where i = Item x (f x)
        nMoves = length $ fromCycle c

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO
