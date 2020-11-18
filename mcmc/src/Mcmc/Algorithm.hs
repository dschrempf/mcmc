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
import Mcmc.Chain.Chain
import Mcmc.Environment
import Mcmc.Monitor
import Mcmc.Proposal
import Numeric.Log

-- | TODO: REFACTOR. Documentation.
class Algorithm a where
  algorithmName :: a -> String

  algorithmIteration :: a -> Int

  -- TODO: Splitmix. Remove IO monad as soon as possible.

  algorithmIterate :: a -> IO a

  algorithmAutoTune :: a -> a

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

  algorithmResetAcceptance :: a -> a

  -- -- | Reset acceptance counts.
  -- mcmcResetA :: Mcmc a ()
  -- mcmcResetA = do
  --   mcmcDebugB "Reset acceptance ratios."
  --   s <- get
  --   let a = acceptance s
  --   put $ s {acceptance = resetA a}

  algorithmSummarizeCycle :: a -> BL.ByteString

  -- -- | Print short summary of 'Proposal's in 'Cycle'. See 'summarizeCycle'.
  -- mcmcSummarizeCycle :: Mcmc a BL.ByteString
  -- mcmcSummarizeCycle = do
  --   a <- gets acceptance
  --   c <- gets cycle
  --   return $ summarizeCycle a c

  algorithmOpenMonitors :: a -> IO ()

  -- -- Monitor.
  -- let m = monitor s
  --     n = iteration s
  --     nm = name s
  -- frc <- reader overwrite
  -- m' <- if n == 0 then liftIO $ mOpen nm frc m else liftIO $ mAppend nm m
  -- put $ s {monitor = m', start = Just (n, t)}

  algorithmExecuteMonitors :: Environment -> a -> IO ()

  -- -- | Execute the 'Monitor's of the chain. See 'mExec'.
  -- mcmcMonitorExec :: ToJSON a => Mcmc a ()
  -- mcmcMonitorExec = do
  --   vb <- reader (verbosity . settings)
  --   s <- get
  --   let i = iteration s
  --       j = iterations s + fromMaybe 0 (burnInIterations s)
  --       m = monitor s
  --       (ss, st) = fromMaybe (error "mcmcMonitorExec: Starting state and time not set.") (start s)
  --       tr = trace s
  --   mt <- liftIO $ mExec vb i ss st tr j m
  --   forM_ mt mcmcOutB

  algorithmCloseMonitors :: a -> IO ()

  -- -- Close the 'Monitor's of the chain. See 'mClose'.
  -- mcmcClose :: ToJSON a => Mcmc a ()
  -- mcmcClose = do
  --   s <- get
  --   mcmcSummarizeCycle >>= mcmcInfoB
  --   mcmcInfoB "Metropolis-Hastings sampler finished."
  --   let m = monitor s
  --   m' <- liftIO $ mClose m
  --   put $ s {monitor = m'}
  --   mcmcSave
  --   t <- liftIO getCurrentTime
  --   let rt = case start s of
  --         Nothing -> error "mcmcClose: Start time not set."
  --         Just (_, st) -> t `diffUTCTime` st
  --   mcmcInfoB $ "Wall clock run time: " <> renderDuration rt <> "."
  --   mcmcInfoS $ "End time: " <> fTime t
  --   case logHandle s of
  --     Just h -> liftIO $ hClose h
  -- Nothing -> return ()

  -- | Save chain(s) with trace of given maximum length.
  algorithmSaveWith :: Int -> a -> BL.ByteString

  -- | Load chain(s).
  algorithmLoadWith ::
    PriorFunction a ->
    LikelihoodFunction a ->
    Cycle a ->
    Monitor a ->
    BL.ByteString ->
    a

  -- | Report prior and likelihood; useful for debugging.
  algorithmReport :: a -> (Log Double, Log Double)
