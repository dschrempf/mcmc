{- |
Module      :  Statistics.Mcmc.Types
Description :  What is an MCMC?
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 18:01:15 2020.

-}

module Statistics.Mcmc.Types
  ( State (..)
  , Mcmc (..)
  , Link (..)
  , Chain (..)
  ) where

import Numeric.Log

-- | A state in the parameter space @a@.
newtype State a = State a
  deriving (Show, Read)

data Mcmc a = Mcmc
  {
    -- ^ The present state.
    state :: State a
    -- ^ The log-likelihood function.
  , llhf  :: State a -> Log Double
    -- We could think of adding a set of moves, or something equivalent, but we
    -- don't do that at the moment :).
  }

-- | A link in the chain.
data Link a = Link
  {
    -- ^ The present state.
    lnState :: State a
    -- ^ The present log-likelihood.
  , lnLlh   :: Log Double
  }
  deriving (Show, Read)

-- | A chain passes through a list of states.
newtype Chain a = Chain [Link a]
  deriving (Show, Read)
