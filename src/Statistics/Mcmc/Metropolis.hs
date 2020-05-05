{- |
Module      :  Statistics.Mcmc.Metropolis
Description :  Metropolis-Hastings at its best
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 20:11:30 2020.

-}

module Statistics.Mcmc.Metropolis
  (
  ) where

import Statistics.Mcmc.Types

mhMove :: Mcmc a -> Move m a -> Link a
mhMove (Mcmc s f) (Move p q) = undefined
  -- 1. Sample new state.
  -- 2. Calculate Metropolis-Hastings ratio.
  -- 3. Accept or reject.
