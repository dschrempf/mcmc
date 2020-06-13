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
  ( McmcData(..)
  , Mcmc
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
import qualified Data.Map.Strict               as M
import           Data.Maybe
import qualified Data.Text.Lazy.IO             as T
import           Data.Time.Clock
import           Data.Time.Format

import           Lens.Micro

import           Mcmc.Item
import           Mcmc.Status
import           Mcmc.Monitor
import           Mcmc.Monitor.Time
import           Mcmc.Move

data McmcData a = McmcData
  {
    mcmcStatus :: Status a
  , mcmcSpec   :: Spec a
  }

lst :: Lens' (McmcData a) (Status a)
lst = lens mcmcStatus (\m s' -> m { mcmcStatus = s' })

lsp :: Lens' (McmcData a) (Spec a)
lsp = lens mcmcSpec (\m s' -> m { mcmcSpec = s' })

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (McmcData a) IO

-- Tune the 'Move's in the 'Cycle' of the Markov chain 'Status'; check
-- acceptance ratio of the last n moves. Tuning has no effect on 'Move's that
-- cannot be tuned. See 'autotune'.
autotuneS :: Int -> Status a -> Spec a -> Spec a
autotuneS n st sp = sp { cycle = mapCycle tuneF (cycle sp) }
 where
  ars = acceptanceRatios n $ acceptance st
  tuneF m = fromMaybe m (autotune (ars M.! m) m)

-- | Auto tune the 'Move's in the 'Cycle' of the chain. See 'autotune'.
mcmcAutotune :: Int -> Mcmc a ()
mcmcAutotune t = do
  d <- get
  let
    st = d^.lst
    sp = d^.lsp
  put $ lsp .~ autotuneS t st sp $ d

-- | Print short summary of 'Move's in 'Cycle'. See 'summarizeCycle'.
mcmcSummarizeCycle :: Maybe Int -> Mcmc a ()
mcmcSummarizeCycle Nothing = do
  liftIO $ putStrLn ""
  c <- gets $ cycle . mcmcSpec
  liftIO $ putStr $ summarizeCycle Nothing c
  liftIO $ putStrLn ""
mcmcSummarizeCycle (Just n) = do
  liftIO $ putStrLn ""
  a <- gets $ acceptance . mcmcStatus
  c <- gets $ cycle . mcmcSpec
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
  d  <- get
  let m = monitor $ d^.lsp
  m' <- liftIO $ mOpen m
  put $ lsp %~ (\sp -> sp { monitor = m', starttime = Just t, totalIterations = Just n }) $ d

-- | Print header line of 'Monitor' (only standard output).
mcmcMonitorHeader :: Mcmc a ()
mcmcMonitorHeader = do
  m <- gets $ monitor . mcmcSpec
  liftIO $ mHeader m

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: Mcmc a ()
mcmcMonitorExec = do
  d  <- get
  ct <- liftIO getCurrentTime
  let sp = mcmcSpec d
      st = mcmcStatus d
      i = iteration st
      j = fromMaybe
        (error "mcmcMonitorExec: Total number of iterations not set.")
        (totalIterations sp)
      (Item x p l) = item st
      m            = monitor sp
      mst          = starttime sp
      dt            = case mst of
        Nothing -> error "mcmcMonitorExec: Start time not set."
        Just t -> ct `diffUTCTime` t
  liftIO $ mExec i p l dt x j m

-- | Close the 'Monitor's of the chain. See 'mClose'.
mcmcClose :: Mcmc a ()
mcmcClose = do
  sp <- gets mcmcSpec
  let m = monitor sp
  liftIO $ mClose m
  t <- liftIO getCurrentTime
  let rt = case starttime sp of
        Nothing -> error "mcmcClose: Start time not set."
        Just st -> t `diffUTCTime` st
  liftIO $ T.putStrLn $ "-- Wall clock run time: " <> renderDuration rt <> "."
  liftIO $ putStrLn $ "-- End time: " <> fTime t
