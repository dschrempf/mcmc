{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Mcmc
-- Description :  Mcmc helpers
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 10:19:45 2020.
module Mcmc.Mcmc
  ( Mcmc,
    mcmcAutotune,
    mcmcSummarizeCycle,
    mcmcInit,
    mcmcReport,
    mcmcMonitorHeader,
    mcmcMonitorExec,
    mcmcClose,
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.State.Strict
import Data.Aeson
import Data.Maybe
import qualified Data.Text.Lazy.IO as T
import Data.Time.Clock
import Data.Time.Format
import Mcmc.Monitor
import Mcmc.Monitor.Time
import Mcmc.Move
import Mcmc.Save
import Mcmc.Status
import Prelude hiding (cycle)

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
  put $ s {cycle = c'}

-- | Print short summary of 'Move's in 'Cycle'. See 'summarizeCycle'.
mcmcSummarizeCycle :: Maybe Int -> Mcmc a ()
mcmcSummarizeCycle Nothing = do
  liftIO $ putStrLn ""
  c <- gets cycle
  liftIO $ T.putStr $ summarizeCycle Nothing c
  liftIO $ putStrLn ""
mcmcSummarizeCycle (Just n) = do
  liftIO $ putStrLn ""
  a <- gets acceptance
  c <- gets cycle
  liftIO $ T.putStr $ summarizeCycle (Just (n, a)) c
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
  let m = monitor s
      n = iteration s
      nm = name s
  m' <- if n == 0
        then liftIO $ mOpen nm m
        else liftIO $ mAppend nm m
  put $ s {monitor = m', start = Just (n, t)}

-- | Report what is going to be done.
mcmcReport :: Mcmc a ()
mcmcReport = do
  s <- get
  let b = burnInIterations s
      t = autoTuningPeriod s
      n = iterations s
  case b of
    Just b' -> liftIO $ putStrLn $ "-- Burn in for " <> show b' <> " iterations."
    Nothing -> return ()
  case t of
    Just t' ->
      liftIO $ putStrLn $
        "-- Auto tune every "
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
  liftIO $ if save s
    then do putStrLn "-- Save Markov chain. For long chains, this may take a while."
            saveStatus (name s <> ".mcmc") s
            putStrLn "-- Done saving Markov chain."
    else putStrLn "-- Do not save the Markov chain."

-- | Execute the 'Monitor's of the chain. See 'mExec'.
mcmcMonitorExec :: ToJSON a => Mcmc a ()
mcmcMonitorExec = do
  s <- get
  let i = iteration s
      j = iterations s + fromMaybe 0 (burnInIterations s)
      m = monitor s
      st = fromMaybe (error "mcmcMonitorExec: Starting state and time not set.") (start s)
      tr = trace s
  liftIO $ mExec i st tr j m

-- | Close the 'Monitor's of the chain. See 'mClose'.
mcmcClose :: ToJSON a => Mcmc a ()
mcmcClose = do
  s <- get
  let n = iterations s
  mcmcMonitorExec
  mcmcSummarizeCycle (Just n)
  liftIO $ putStrLn "-- Metropolis-Hastings sampler finished."
  let m = monitor s
  m' <- liftIO $ mClose m
  put $ s {monitor = m'}
  mcmcSave
  t <- liftIO getCurrentTime
  let rt = case start s of
        Nothing -> error "mcmcClose: Start time not set."
        Just (_, st) -> t `diffUTCTime` st
  liftIO $ T.putStrLn $ "-- Wall clock run time: " <> renderDuration rt <> "."
  liftIO $ putStrLn $ "-- End time: " <> fTime t
