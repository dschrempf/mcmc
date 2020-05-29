{- |
Module      :  Statistics.Mcmc.Mcmc
Description :  Mcmc helpers
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May 29 10:19:45 2020.

-}

module Statistics.Mcmc.Mcmc
  (
    Mcmc
  , mcmcAutotune
  , mcmcSummarizeCycle
  , mcmcMonitorOpen
  , mcmcMonitorHeader
  , mcmcMonitorExec
  , mcmcMonitorClose
  ) where

import Prelude hiding (cycle)

import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict hiding (state)
import qualified Data.Map.Strict as M
import Data.Maybe

import Statistics.Mcmc.Item
import Statistics.Mcmc.Status
import Statistics.Mcmc.Monitor
import Statistics.Mcmc.Move

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO

-- Tune the 'Move's in the 'Cycle' of the Markov chain 'Status'; check
-- acceptance ratio of the last n moves. Tuning has no effect on 'Move's that
-- cannot be tuned. See 'autotune'.
autotuneS :: Int -> Status a -> Status a
autotuneS n s = s {cycle = mapCycle tuneF (cycle s)}
  where
    ars  = acceptanceRatios n $ acceptance s
    tuneF m = fromMaybe m (autotune (ars M.! m) m)

-- | Auto tune the 'Move's in the 'Cycle' of the chain. See 'autotune'.
mcmcAutotune :: Int -> Mcmc a ()
mcmcAutotune t = modify' (autotuneS t)

-- | Print short summary of 'Move's in 'Cycle'. See 'summarizeCycle'.
mcmcSummarizeCycle :: Maybe Int -> Mcmc a ()
mcmcSummarizeCycle Nothing = do
  liftIO $ putStrLn ""
  c <- gets cycle
  liftIO $ putStr $ summarizeCycle Nothing c
  liftIO $ putStrLn ""
mcmcSummarizeCycle (Just n) = do
  liftIO $ putStrLn ""
  a <- gets acceptance
  c <- gets cycle
  liftIO $ putStr $ summarizeCycle (Just (n, a)) c
  liftIO $ putStrLn ""

-- | Open the 'Monitor's of the chain. See 'mOpen'.
mcmcMonitorOpen :: Mcmc a ()
mcmcMonitorOpen = do
  s  <- get
  m  <- gets monitor
  m' <- liftIO $ mOpen m
  put s { monitor = m' }

-- | Print header line of 'Monitor' (only standard output).
mcmcMonitorHeader :: Mcmc a ()
mcmcMonitorHeader = do
  m <- gets monitor
  liftIO $ mHeader m

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: Mcmc a ()
mcmcMonitorExec = do
  s <- get
  let i = iteration s
      x = state $ item s
      m = monitor s
  liftIO $ mExec i x m

-- | Close the 'Monitor's of the chain. See 'mClose'.
mcmcMonitorClose :: Mcmc a ()
mcmcMonitorClose = do
  m <- gets monitor
  liftIO $ mClose m
