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
--
-- Functions to work with the 'Mcmc' state transformer.
module Mcmc.Mcmc
  ( mcmc,
  )
where

import Codec.Compression.GZip
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.RWS.CPS
import Data.Aeson
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Time.Clock
import Data.Time.Format
import Mcmc.Algorithm
import Mcmc.Environment
import Mcmc.Monitor.Time
import Mcmc.Settings
import System.IO
import Prelude hiding (cycle)

-- The MCMC algorithm has read access to an environment and uses an algorithm
-- transforming the state @a@.
type Mcmc a = RWST Environment () a IO

msgPrepare :: Char -> BL.ByteString -> BL.ByteString
msgPrepare c t = BL.cons c $ ": " <> t

-- Write to standard output and log file.
mcmcOutB :: BL.ByteString -> Mcmc a ()
mcmcOutB msg = do
  h <- reader logHandle
  liftIO $ BL.putStrLn msg >> BL.hPutStrLn h msg

-- -- Perform warning action.
-- mcmcWarnA :: Mcmc a () -> Mcmc a ()
-- mcmcWarnA a = reader (verbosity . settings) >>= \v -> when (v >= Warn) a

-- -- Print warning message.
-- mcmcWarnB :: BL.ByteString -> Mcmc a ()
-- mcmcWarnB = mcmcWarnA . mcmcOutB . msgPrepare 'W'

-- -- Print warning message.
-- mcmcWarnS :: String -> Mcmc a ()
-- mcmcWarnS = mcmcWarnB . BL.pack

-- Perform info action.
mcmcInfoA :: Mcmc a () -> Mcmc a ()
mcmcInfoA a = reader (verbosity . settings) >>= \v -> when (v >= Info) a

-- Print info message.
mcmcInfoB :: BL.ByteString -> Mcmc a ()
mcmcInfoB = mcmcInfoA . mcmcOutB . msgPrepare 'I'

-- Print info message.
mcmcInfoS :: String -> Mcmc a ()
mcmcInfoS = mcmcInfoB . BL.pack

-- Perform debug action.
mcmcDebugA :: Mcmc a () -> Mcmc a ()
mcmcDebugA a = reader (verbosity . settings) >>= \v -> when (v == Debug) a

-- Print debug message.
mcmcDebugB :: BL.ByteString -> Mcmc a ()
mcmcDebugB = mcmcDebugA . mcmcOutB . msgPrepare 'D'

-- Print debug message.
mcmcDebugS :: String -> Mcmc a ()
mcmcDebugS = mcmcDebugB . BL.pack

fTime :: FormatTime t => t -> String
fTime = formatTime defaultTimeLocale "%B %-e, %Y, at %H:%M %P, %Z."

mcmcReportTime :: Mcmc a ()
mcmcReportTime = do
  ti <- reader startingTime
  mcmcInfoS $ "Starting time: " <> fTime ti

mcmcExecute :: Algorithm a => Mcmc a ()
mcmcExecute = do
  s <- reader settings
  let iTotal = iterations s
  when (iTotal <= 0) $
    error "mcmcExecute: Total number of iterations is zero or negative."
  case executionMode s of
    Fail -> mcmcNewRun
    Overwrite -> mcmcNewRun
    Continue -> mcmcContinue

mcmcNewRun :: Algorithm a => Mcmc a ()
mcmcNewRun = do
  s <- reader settings
  mcmcInfoB "Start new MCMC sampler."
  mcmcInfoB "Initial state."
  mcmcExecuteMonitors
  get >>= mcmcInfoB . algorithmSummarizeCycle
  mcmcBurnIn
  mcmcResetAcceptance
  let i = iterations s
  mcmcInfoS $ "Run chain for " ++ show i ++ " iterations."
  mcmcIterate i

mcmcContinue :: Algorithm a => Mcmc a ()
mcmcContinue = do
  s <- reader settings
  let iTotal = iterations s
  mcmcInfoB "Continuation of MCMC sampler."
  a <- get
  let iCurrent = algorithmIteration a
  mcmcInfoS $ "Current iteration: " ++ show iCurrent ++ "."
  mcmcInfoS $ "Total iterations: " ++ show iTotal ++ "."
  let di = iTotal - iCurrent
  when (di <= 0) $
    error "mcmcContinue: Current iteration is equal or larger to the total number of iterations."
  get >>= mcmcInfoB . algorithmSummarizeCycle
  mcmcInfoS $ "Run chain for " ++ show di ++ " iterations."
  mcmcIterate di

mcmcBurnIn :: Algorithm a => Mcmc a ()
mcmcBurnIn = do
  s <- reader settings
  case burnIn s of
    NoBurnIn ->
      mcmcInfoS "No burn in."
    BurnInNoAutoTuning n -> do
      mcmcInfoS $ "Burn in for " <> show n <> " iterations."
      when (n < 0) $ error "mcmcBurnIn: Number of burn in iterations is negative."
      mcmcInfoS "Auto tuning is disabled."
      mcmcIterate n
      get >>= mcmcInfoB . algorithmSummarizeCycle
      mcmcInfoB "Burn in finished."
    BurnInWithAutoTuning n t -> do
      mcmcInfoS $ "Burn in for " ++ show n ++ " iterations."
      when (n < 0) $ error "mcmcBurnIn: Number of burn in iterations is negative."
      mcmcInfoS $ "Auto tuning is enabled with a period of "++ show t ++ "."
      when (t <= 0) $ error "mcmcBurnIn: Auto tuning period is zero or negative."
      mcmcBurnInWithAutoTuning n t
      mcmcInfoB "Burn in finished."

mcmcBurnInWithAutoTuning :: Algorithm a => Int -> Int -> Mcmc a ()
mcmcBurnInWithAutoTuning b t
  | b > t = do
    mcmcResetAcceptance
    mcmcIterate t
    get >>= mcmcDebugB . algorithmSummarizeCycle
    mcmcAutotune
    mcmcBurnInWithAutoTuning (b - t) t
  | otherwise = do
    mcmcResetAcceptance
    mcmcIterate b
    get >>= mcmcInfoB . algorithmSummarizeCycle
    mcmcInfoS $ "Acceptance ratios calculated over the last " <> show b <> " iterations."

mcmcIterate :: Algorithm a => Int -> Mcmc a ()
mcmcIterate n | n < 0 = error "mcmcIterate: Number of iterations is negative."
              | n == 0 = return ()
              | otherwise = do
                  -- TODO: Splitmix. Remove IO monad as soon as possible.
                  a' <- get >>= liftIO . algorithmIterate
                  put a'
                  mcmcExecuteMonitors
                  mcmcIterate (n-1)

-- Auto tune the proposals.
mcmcAutotune :: Algorithm a => Mcmc a ()
mcmcAutotune = do
  mcmcDebugB "Auto tune."
  modify algorithmAutoTune

-- Reset acceptance counts.
mcmcResetAcceptance :: Algorithm a => Mcmc a ()
mcmcResetAcceptance = do
  mcmcDebugB "Reset acceptance ratios."
  modify algorithmResetAcceptance

-- Execute the monitors of the chain.
mcmcExecuteMonitors :: Algorithm a => Mcmc a ()
mcmcExecuteMonitors = do
  e <- ask
  a <- get
  liftIO $ algorithmExecuteMonitors e a

-- Save the MCMC run.
mcmcSave :: Algorithm a => Mcmc a ()
mcmcSave = do
  ss <- reader settings
  a <- get
  case saveMode ss of
    NoSave -> mcmcInfoB "Do not save the Markov chain."
    SaveWithTrace n -> do
      let fns = name ss ++ ".settings"
      mcmcInfoS $ "Save settings to " ++ fns ++ "."
      liftIO $ BL.writeFile fns $ encode ss
      let fnc = name ss ++ ".chain"
      mcmcInfoS $
        "Save compressed Markov chain with trace of length "
          ++ show n
          ++ " to "
          ++ fnc
          ++ "."
      mcmcInfoB "For long traces, or complex objects, this may take a while."
      liftIO $ BL.writeFile fnc $ compress $ algorithmSaveWith n a
      mcmcInfoB "Markov chain saved."

-- Report and finish up.
mcmcClose :: Algorithm a => Mcmc a ()
mcmcClose = do
  get >>= mcmcInfoB . algorithmSummarizeCycle
  mcmcInfoB "Metropolis-Hastings sampler finished."
  mcmcSave
  ti <- reader startingTime
  te <- liftIO getCurrentTime
  let dt = te `diffUTCTime` ti
  mcmcInfoB $ "Wall clock run time: " <> renderDuration dt <> "."
  mcmcInfoS $ "End time: " <> fTime te
  h <- reader logHandle
  liftIO $ hClose h

-- Initialize the run, execute the run, and close the run.
mcmcRun :: Algorithm a => Mcmc a ()
mcmcRun = do
  mcmcDebugB "The settings are:"
  reader settings >>= mcmcDebugS . show

  -- Initialize.
  get >>= mcmcInfoS . algorithmName
  get >>= liftIO . algorithmOpenMonitors
  mcmcReportTime

  -- Execute.
  mcmcExecute

  -- Close.
  mcmcClose
  get >>= liftIO . algorithmCloseMonitors

-- | Run an MCMC algorithm.
mcmc :: Algorithm a => Settings -> a -> IO a
mcmc s a = do
  e <- initializeEnvironment s
  fst <$> execRWST mcmcRun e a
