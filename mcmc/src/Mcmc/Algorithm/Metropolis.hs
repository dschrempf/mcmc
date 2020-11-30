{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Algorithm.Metropolis
-- Description :  Metropolis-Hastings-Green algorithm
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 20:11:30 2020.
module Mcmc.Algorithm.Metropolis
  ( MHG (..),
    mhg,
    mhgSave,
    mhgLoad,
    mhgAccept,
  )
where

import Codec.Compression.GZip
import Control.Monad
import Control.Monad.IO.Class
import Data.Aeson
import qualified Data.ByteString.Lazy.Char8 as BL
import Mcmc.Algorithm
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Save
import Mcmc.Chain.Trace
import Mcmc.Environment
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log
import System.Random.MWC
import Prelude hiding (cycle)

-- | The Metropolis-Hastings-Green (MHG) algorithm.
newtype MHG a = MHG {fromMHG :: Chain a}

instance ToJSON a => Algorithm (MHG a) where
  aName = const "Metropolis-Hastings-Green (MHG)"
  aIteration = iteration . fromMHG
  aIterate = mhgIterate
  aAutoTune = mhgAutoTune
  aResetAcceptance = mhgResetAcceptance
  aSummarizeCycle = mhgSummarizeCycle
  aOpenMonitors = mhgOpenMonitors
  aExecuteMonitors = mhgExecuteMonitors
  aStdMonitorHeader = mhgStdMonitorHeader
  aCloseMonitors = mhgCloseMonitors
  aSave = mhgSave

-- | Initialize an MHG algorithm.
mhg ::
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  -- | The initial state in the state space @a@.
  a ->
  -- | A source of randomness. For reproducible runs, make sure to use
  -- generators with the same, fixed seed.
  GenIO ->
  MHG a
mhg pr lh cc mn i0 g = MHG $ Chain l0 0 tr ac g 0 pr lh cc mn
  where
    l0 = Link i0 (pr i0) (lh i0)
    tr = singletonT l0
    ac = emptyA $ ccProposals cc

mhgFn :: String -> FilePath
mhgFn nm = nm ++ ".chain"

-- | Save an MHG algorithm.
mhgSave ::
  ToJSON a =>
  -- | Maximum length of trace.
  Int ->
  -- | Analysis name.
  String ->
  MHG a ->
  IO ()
mhgSave n nm (MHG c) = BL.writeFile (mhgFn nm) $ compress $ encode $ toSavedChain n c

-- | Load an MHG algorithm.
mhgLoad ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  -- | Analysis name.
  String ->
  IO (MHG a)
mhgLoad pr lh cc mn nm =
  MHG
    . either error (fromSavedChain pr lh cc mn)
    . eitherDecode
    . decompress
    <$> BL.readFile (mhgFn nm)

-- The MHG ratio.
--
-- 'Infinity' if fX is zero. In this case, the proposal is always accepted.
--
-- 'NaN' if (fY or q) and fX are zero. In this case, the proposal is always
-- rejected.

-- There is a discrepancy between authors saying that one should (a) always
-- accept the new state when the current posterior is zero (Chapter 4 of the
-- Handbook of Markov Chain Monte Carlo), or (b) almost surely reject the
-- proposal when either fY or q are zero (Chapter 1). Since I trust the author
-- of Chapter 1 (Charles Geyer) I choose to follow option (b).
mhgRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
-- q = qYX / qXY * jXY; see 'ProposalSimple'.
-- j = Jacobian.
mhgRatio fX fY q j = fY / fX * q * j
{-# INLINE mhgRatio #-}

-- | Accept or reject a proposal with given MHG ratio?
mhgAccept :: Log Double -> GenIO -> IO Bool
mhgAccept r g
  | ln r >= 0.0 = return True
  | otherwise = do
    b <- uniform g
    return $ b < exp (ln r)

mhgPropose :: MHG a -> Proposal a -> IO (MHG a)
mhgPropose (MHG c) p = do
  -- 1. Sample new state.
  (!y, !q, !j) <- liftIO $ s x g
  -- 2. Calculate Metropolis-Hastings-Green ratio.
  let !pY = pF y
      !lY = lF y
      !r = mhgRatio (pX * lX) (pY * lY) q j
  -- 3. Accept or reject.
  -- if ln r >= 0.0
  --   then do
  --     let !ac' = pushA p True ac
  --     return $ MHG $ c {link = Link y pY lY, acceptance = ac'}
  --   else do
  --     b <- uniform g
  --     if b < exp (ln r)
  --       then do
  --         let !ac' = pushA p True ac
  --         return $ MHG $ c {link = Link y pY lY, acceptance = ac'}
  --       else do
  --         let !ac' = pushA p False ac
  --         return $ MHG $ c {acceptance = pushA p False ac'}
  accept <- mhgAccept r g
  if accept
    then do
      let !ac' = pushA p True ac
      return $ MHG $ c {link = Link y pY lY, acceptance = ac'}
    else do
      let !ac' = pushA p False ac
      return $ MHG $ c {acceptance = pushA p False ac'}
  where
    s = pSimple p
    (Link x pX lX) = link c
    pF = priorFunction c
    lF = likelihoodFunction c
    ac = acceptance c
    g = generator c

mhgPush :: MHG a -> MHG a
mhgPush (MHG c) = MHG c {trace = pushT i t, iteration = succ n}
  where
    i = link c
    t = trace c
    n = iteration c

-- Ignore the number of capabilities. I have tried a lot of stuff, but the MHG
-- algorithm is just inherently sequential. Parallelization can be achieved by
-- having parallel prior and/or likelihood functions, or by using algorithms
-- running parallel chains such as 'MC3'.
mhgIterate :: Int -> MHG a -> IO (MHG a)
mhgIterate _ a = do
  ps <- orderProposals cc g
  a' <- foldM mhgPropose a ps
  return $ mhgPush a'
  where
    c = fromMHG a
    cc = cycle c
    g = generator c

mhgAutoTune :: MHG a -> MHG a
mhgAutoTune (MHG c) = MHG $ c {cycle = autoTuneCycle ac cc}
  where
    ac = acceptance c
    cc = cycle c

mhgResetAcceptance :: MHG a -> MHG a
mhgResetAcceptance (MHG c) = MHG $ c {acceptance = resetA ac}
  where
    ac = acceptance c

mhgSummarizeCycle :: MHG a -> BL.ByteString
mhgSummarizeCycle (MHG c) = summarizeCycle ac cc
  where
    cc = cycle c
    ac = acceptance c

mhgOpenMonitors :: Environment -> MHG a -> IO (MHG a)
mhgOpenMonitors e (MHG c) = do
  m' <- mOpen nm em m
  return $ MHG c {monitor = m'}
  where
    m = monitor c
    s = settings e
    nm = sAnalysisName s
    em = sExecutionMode s

mhgExecuteMonitors :: Environment -> MHG a -> IO (Maybe BL.ByteString)
mhgExecuteMonitors e (MHG c) = mExec vb i i0 t0 tr j m
  where
    s = settings e
    vb = sVerbosity s
    i = iteration c
    i0 = start c
    t0 = startingTime e
    tr = trace c
    b = burnInIterations $ sBurnIn s
    j = sIterations s + b
    m = monitor c

mhgStdMonitorHeader :: MHG a -> BL.ByteString
mhgStdMonitorHeader (MHG c) = msHeader (mStdOut $ monitor c)

mhgCloseMonitors :: MHG a -> IO (MHG a)
mhgCloseMonitors (MHG c) = do
  m' <- mClose m
  return $ MHG $ c {monitor = m'}
  where
    m = monitor c
