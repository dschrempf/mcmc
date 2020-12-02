-- |
-- Module      :  Mcmc.Chain.Chain
-- Description :  Simple representation of a Markov chain
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 18:01:15 2020.
module Mcmc.Chain.Chain
  ( PriorFunction,
    noPrior,
    LikelihoodFunction,
    noLikelihood,
    Chain (..),
  )
where

-- Note: It is not necessary to add another type @b@ to store supplementary
-- information about the chain. The information can just be stored in @a@
-- equally well.

import Mcmc.Chain.Link
import Mcmc.Chain.Trace
import Mcmc.Monitor
import Mcmc.Proposal
import Numeric.Log
import System.Random.MWC hiding (save)
import Prelude hiding (cycle)

-- | Prior function.
type PriorFunction a = a -> Log Double

-- | Flat prior function. Useful for testing and debugging.
noPrior :: PriorFunction a
noPrior = const 1.0

-- | Likelihood function.
type LikelihoodFunction a = a -> Log Double

-- | Flat likelihood function. Useful for testing and debugging.
noLikelihood :: LikelihoodFunction a
noLikelihood = const 1.0

-- | The chain contains all information to run an MCMC sampler.
--
-- The state of a chain has type @a@. If necessary, the type @a@ can also be
-- used to store auxiliary information.
--
-- For example, the chain stores information about the current 'Link' and
-- 'iteration', the 'Trace', the 'Acceptance' rates, and the random number
-- generator.
--
-- Further, the chain includes auxiliary variables and functions such as the
-- prior and likelihood functions, or 'Proposal's to move around the state space
-- and to 'Monitor' an MCMC run.
--
-- The 'Mcmc.Environment.Environment' of the chain is not stored externally.
data Chain a = Chain
  { -- Variables; saved.
    -- | Chain index; useful if more chains are run.
    chainId :: Int,
    -- | The current 'Link' of the chain combines the current state and the
    -- current likelihood. The link is updated after a proposal has been
    -- executed.
    link :: Link a,
    -- | The current iteration or completed number of cycles.
    iteration :: Int,
    -- | The 'Trace' of the Markov chain in reverse order, the most recent
    -- 'Link' is at the head of the list. In contrast to the link, the trace is
    -- updated only after all proposals in the cycle have been executed.
    trace :: Trace a,
    -- | For each 'Proposal', store the list of accepted (True) and rejected (False)
    -- proposals; for reasons of efficiency, the list is also stored in reverse
    -- order.
    acceptance :: Acceptance (Proposal a),
    -- | The random number generator.
    generator :: GenIO,
    --
    -- Variables and functions; not saved.

    -- | Starting iteration of the chain; used to calculate run time and ETA.
    start :: Int,
    -- | The prior function. The un-normalized posterior is the product of the
    -- prior and the likelihood.
    priorFunction :: PriorFunction a,
    -- | The likelihood function. The un-normalized posterior is the product of
    -- the prior and the likelihood.
    likelihoodFunction :: LikelihoodFunction a,
    -- | A set of 'Proposal's form a 'Cycle'.
    cycle :: Cycle a,
    -- | A 'Monitor' observing the chain.
    monitor :: Monitor a
  }
