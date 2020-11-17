-- |
-- Module      :  Mcmc.Algorithm
-- Description :  MCMC algorithms
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 16 14:37:11 2020.
module Mcmc.Algorithm
  ( Algorithm (..),
  )
where

import qualified Data.ByteString.Lazy.Char8 as BL
import Mcmc.Chain
import Mcmc.Environment
import Mcmc.Monitor
import Mcmc.Proposal
import Numeric.Log

-- | TODO: REFACTOR. Documentation.
class Algorithm a where
  -- | Get current iteration.
  getIteration :: a -> Int

  -- TODO: Splitmix. Remove IO monad as soon as possible.

  -- | Perform one iteration.
  jump :: a -> IO a

  -- | Auto tune the proposals.
  autoTune :: a -> a

  -- -- | Auto tune the 'Proposal's in the 'Cycle' of the chain. Reset acceptance counts.
  -- -- See 'autoTuneCycle'.
  -- mcmcAutoTune :: Mcmc a ()
  -- mcmcAutoTune = do
  --   mcmcDebugB "Auto tune."
  --   s <- get
  --   let a = acceptance s
  --       c = cycle s
  --       c' = autoTuneCycle a c
  --   put $ s {cycle = c'}

  -- | Reset acceptance ratios.
  resetAcceptance :: a -> a

  -- | Summarize cycle.
  summarizeCycle :: a -> BL.ByteString

  -- | Open the monitor files.
  openMonitors :: a -> IO ()

  -- -- Monitor.
  -- let m = monitor s
  --     n = iteration s
  --     nm = name s
  -- frc <- reader overwrite
  -- m' <- if n == 0 then liftIO $ mOpen nm frc m else liftIO $ mAppend nm m
  -- put $ s {monitor = m', start = Just (n, t)}

  -- | Execute the monitors.
  execMonitors :: Environment -> a -> IO ()

  -- -- | Print short summary of 'Proposal's in 'Cycle'. See 'summarizeCycle'.
  -- mcmcSummarizeCycle :: Mcmc a BL.ByteString
  -- mcmcSummarizeCycle = do
  --   a <- gets acceptance
  --   c <- gets cycle
  --   return $ summarizeCycle a c

  -- -- | Reset acceptance counts.
  -- mcmcResetA :: Mcmc a ()
  -- mcmcResetA = do
  --   mcmcDebugB "Reset acceptance ratios."
  --   s <- get
  --   let a = acceptance s
  --   put $ s {acceptance = resetA a}

  -- | Clean the state.
  clean :: a -> a

  -- -- | Clean the state.
  -- mcmcClean :: Mcmc a ()
  -- mcmcClean = do
  --   s <- get
  --   let cl = cleaner s
  --       i = iteration s
  --   case cl of
  --     Just (Cleaner n f) | i `mod` n == 0 -> do
  --       mcmcDebugB "Clean state."
  --       let (Item st pr lh) = item s
  --       mcmcDebugS $
  --         "Old log prior and log likelihood: " ++ show (ln pr) ++ ", " ++ show (ln lh) ++ "."
  --       let prF = priorF s
  --           lhF = likelihoodF s
  --           st' = f st
  --           pr' = prF st'
  --           lh' = lhF st'
  --       mcmcDebugS $
  --         "New log prior and log likelihood: " ++ show (ln pr') ++ ", " ++ show (ln lh') ++ "."
  --       let dLogPr = abs $ ln pr - ln pr'
  --           dLogLh = abs $ ln lh - ln lh'
  --       when
  --         (dLogPr > 0.01)
  --         (mcmcWarnS $ "Log of old and new prior differ by " ++ show dLogPr ++ ".")
  --       when
  --         (dLogPr > 0.01)
  --         (mcmcWarnS $ "Log of old and new likelihood differ by " ++ show dLogLh ++ ".")
  --       put $ s {item = Item st' pr' lh'}
  --     _ -> return ()

  -- | Save chain(s) with trace of given maximum length.
  saveWith :: Int -> a -> BL.ByteString

  -- | Load chain(s).
  loadWith ::
    PriorFunction a ->
    LikelihoodFunction a ->
    -- CleaningFunction a ->
    Cycle a ->
    Monitor a ->
    BL.ByteString ->
    a

  -- | Report prior and likelihood; useful for debugging.
  report :: a -> (Log Double, Log Double)
