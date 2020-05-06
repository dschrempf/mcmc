{- |
Module      :  Statistics.Mcmc.Metropolis
Description :  Metropolis-Hastings at its best
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 20:11:30 2020.

TODO: Reexport types?

-}

module Statistics.Mcmc.Metropolis
  ( mh
  ) where

import Prelude hiding (cycle)

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict
import Numeric.Log
import System.Random.MWC

import Statistics.Mcmc.Types

-- import Debug.Trace

{-# INLINE mhRatio #-}
mhRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
mhRatio lX lY qXY qYX = lY * qYX / lX / qXY

mhMove :: Show a => Move a -> Mcmc a Bool
mhMove (Move _ p q) = do
  s <- get
  let (Item x lX) = item s
      f           = logPosteriorF s
      g           = generator s
  -- 1. Sample new state.
  y <- liftIO $ p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let
    lY = f y
    r  = mhRatio lX lY (q x y) (q y x)
  -- 3. Accept or reject.
  -- XXX: Can this be improved in terms of speed?
  -- if traceShow (lX, lY, q x y, q y x, r) $ r >= 1.0
  if ln r >= 0.0
    then put s{item = Item y lY} >> return True
    else do b <- uniform g
            -- Only update the 'Item' after a full cycle.
            if b < exp (ln r)
            then put s{item = Item y lY} >> return True
            else return False

mhCycle :: Show a => (Mcmc a) ()
mhCycle = do
  (Cycle mvs) <- moves <$> get
  a <- mapM mhMove mvs
  s <- get
  let i = item s
      t = trace s
      as = acceptance s
      s' = s {trace = prependT i t, acceptance = prependA a as}
  put s'

-- Run a given number of Metropolis-Hastings cycles.
mhRun :: Show a => Int -> Mcmc a ()
mhRun n = replicateM_ n mhCycle

-- | Run a given number of Metropolis-Hastings cycles.
--
-- TODO: Improve documentation.
--
-- TODO: Simplify type for users.
mh :: Show a => Int -> Status a -> IO (Status a)
mh n = execStateT (mhRun n)
