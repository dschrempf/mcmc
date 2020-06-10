{- |
Module      :  Mcmc.Item
Description :  Links of Markov chains
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May 20 09:10:27 2020.

-}

module Mcmc.Item
  ( Item(..)
  )
where

import           Numeric.Log

-- | An 'Item', or link of the Markov chain. For reasons of computational
-- efficiency, each state is associated with the corresponding log-prior and
-- log-likelihood.
data Item a = Item
  {
    -- | The current state in the state space @a@.
    state         :: a
    -- | The current log-prior.
  , logPrior      :: Log Double
    -- | The current log-likelihood.
  , logLikelihood :: Log Double
  }
  deriving (Eq, Ord, Show, Read)