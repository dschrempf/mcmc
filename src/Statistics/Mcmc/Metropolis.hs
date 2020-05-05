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

mhMove :: Move m a -> Mcmc a -> Link a
mhMove (Move p q) (Mcmc s f) = undefined
  -- 1. Sample new state.
  -- 2. Calculate Metropolis-Hastings ratio.
  -- 3. Accept or reject.

mhCycle :: Cycle m a -> Mcmc a -> Link a
-- mhCycle = fold moves in cycle
mhCycle = undefined

mh :: Int -> Cycle m a -> Mcmc a -> Chain a
-- mh n = repeat mhCycle n times
mh = undefined
