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
  , mcmcReport
  , mcmcMonitorHeader
  , mcmcMonitorExec
  , mcmcClose
  )
where

import           Prelude                 hiding ( cycle )

import           Control.Monad
import           Control.Monad.IO.Class
import           Control.Monad.Trans.State.Strict
import           Data.Aeson
import           Data.Maybe
import qualified Data.Text.Lazy.IO             as T
import           Data.Time.Clock
import           Data.Time.Format

import           Mcmc.Monitor
import           Mcmc.Monitor.Time
import           Mcmc.Move
import           Mcmc.Save
import           Mcmc.Status

-- | An Mcmc state transformer; usually fiddling around with this type is not
-- required, but it is used by the different inference algorithms.
type Mcmc a = StateT (Status a) IO

-- | Auto tune the 'Move's in the 'Cycle' of the chain. See 'autotune'.
mcmcAutotune :: Int -> Mcmc a ()
mcmcAutotune t = do
  s <- get
  let a = acceptance s
      c = cycle s
      c' = autotuneC t a c
  put $ s { cycle = c' }

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
mcmcInit :: Mcmc a ()
mcmcInit = do
  t <- liftIO getCurrentTime
  liftIO $ putStrLn $ "-- Start time: " <> fTime t
  s <- get
  let m  = monitor s
      nm = name s
  m' <- liftIO $ mOpen nm m
  put $ s { monitor = m', starttime = Just t}

-- TODO: This is now tuned to MH, which should not be the case.
-- | Report what is going to be done.
mcmcReport :: Mcmc a ()
mcmcReport = do
  s <- get
  let b = burnInIterations s
      t = autoTuningPeriod s
      n = iterations s
  liftIO $ putStrLn "-- Start of Metropolis-Hastings sampler."
  case b of
    Just b' -> liftIO $ putStrLn $ "-- Burn in for " <> show b' <> " iterations."
    Nothing -> return ()
  case t of
    Just t' ->
      liftIO $ putStrLn
        $  "-- Auto tune every "
        <> show t'
        <> " iterations (during burn in only)."
    Nothing -> return ()
  liftIO $ putStrLn $ "-- Run chain for " <> show n <> " iterations."

-- | Print header line of 'Monitor' (only standard output).
mcmcMonitorHeader :: Mcmc a ()
mcmcMonitorHeader = do
  m <- gets monitor
  liftIO $ mHeader m

-- Save the status of an MCMC run. See 'saveStatus'.
mcmcSave :: ToJSON a => Mcmc a ()
mcmcSave = do
  s <- get
  liftIO $ saveStatus (name s <> ".mcmc") s

-- Save state every 10 seconds.
dtSave :: NominalDiffTime
dtSave = 10

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: ToJSON a => Mcmc a ()
mcmcMonitorExec = do
  s  <- get
  ct <- liftIO getCurrentTime
  let i = iteration s
      j = iterations s + fromMaybe 0 (burnInIterations s)
      m   = monitor s
      mst = starttime s
      dt  = case mst of
        Nothing -> error "mcmcMonitorExec: Start time not set."
        Just st -> ct `diffUTCTime` st
      tr = trace s
  liftIO $ mExec i dt tr j m
  -- TODO: The save should not be here, but at the moment it is convenient.
  case savetime s of
    Nothing  -> mcmcSave >> put s { savetime = Just ct }
    Just svt -> when (ct `diffUTCTime` svt > dtSave)
                     (mcmcSave >> put s { savetime = Just ct })


-- | Close the 'Monitor's of the chain. See 'mClose'.
mcmcClose :: Mcmc a ()
mcmcClose = do
  s <- get
  let m = monitor s
  m' <- liftIO $ mClose m
  put $ s { monitor = m' }
  t <- liftIO getCurrentTime
  let rt = case starttime s of
        Nothing -> error "mcmcClose: Start time not set."
        Just st -> t `diffUTCTime` st
  liftIO $ T.putStrLn $ "-- Wall clock run time: " <> renderDuration rt <> "."
  liftIO $ putStrLn $ "-- End time: " <> fTime t
