{-# LANGUAGE OverloadedStrings #-}

{- |
Module      :  Mcmc.Mcmc
Description :  Mcmc helpers
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May 29 10:19:45 2020.

-}

module Mcmc.Mcmc
  ( Mcmc
  , mcmcAutotune
  , mcmcSummarizeCycle
  , mcmcInit
  , mcmcMonitorHeader
  , mcmcMonitorExec
  , mcmcClose
  )
where

import           Prelude                 hiding ( cycle )

import           Control.Monad.IO.Class
import           Control.Monad.Trans.State.Strict
                                         hiding ( state )
import qualified Data.Map.Strict               as M
import           Data.Maybe
import qualified Data.Text.Lazy.IO             as T
import           Data.Time.Clock
import           Data.Time.Format

import           Mcmc.Item
import           Mcmc.Status
import           Mcmc.Monitor
import           Mcmc.Monitor.Time
import           Mcmc.Move

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO

-- Tune the 'Move's in the 'Cycle' of the Markov chain 'Status'; check
-- acceptance ratio of the last n moves. Tuning has no effect on 'Move's that
-- cannot be tuned. See 'autotune'.
autotuneS :: Int -> Status a -> Status a
autotuneS n s = s { cycle = mapCycle tuneF (cycle s) }
 where
  ars = acceptanceRatios n $ acceptance s
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

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

-- | Set the total number of iterations, the current time and open the
-- 'Monitor's of the chain. See 'mOpen'.
mcmcInit :: Int -> Mcmc a ()
mcmcInit n = do
  t <- liftIO getCurrentTime
  liftIO $ putStrLn $ "-- Start time: " <> fTime t
  s  <- get
  m  <- gets monitor
  m' <- liftIO $ mOpen m
  put s { monitor = m', starttime = Just t, totalIterations = Just n }

-- | Print header line of 'Monitor' (only standard output).
mcmcMonitorHeader :: Mcmc a ()
mcmcMonitorHeader = do
  m <- gets monitor
  liftIO $ mHeader m

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: Mcmc a ()
mcmcMonitorExec = do
  s  <- get
  ct <- liftIO getCurrentTime
  let i = iteration s
      j = fromMaybe
        (error "mcmcMonitorExec: Total number of iterations not set.")
        (totalIterations s)
      (Item x p l) = item s
      m            = monitor s
      mst          = starttime s
      t            = case mst of
        Nothing -> error "mcmcMonitorExec: Start time not set."
        Just st -> ct `diffUTCTime` st
  liftIO $ mExec i p l t x j m

-- | Close the 'Monitor's of the chain. See 'mClose'.
mcmcClose :: Mcmc a ()
mcmcClose = do
  s <- get
  let m = monitor s
  liftIO $ mClose m
  t <- liftIO getCurrentTime
  let rt = case starttime s of
        Nothing -> error "mcmcClose: Start time not set."
        Just st -> t `diffUTCTime` st
  liftIO $ T.putStrLn $ "-- Wall clock run time: " <> renderDuration rt <> "."
  liftIO $ putStrLn $ "-- End time: " <> fTime t
