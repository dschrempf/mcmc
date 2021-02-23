{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Algorithm.MHG
-- Description :  Metropolis-Hastings-Green algorithm
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 20:11:30 2020.
--
-- The Metropolis-Hastings-Green ('MHG') algorithm.
--
-- For example, see Geyer, C. J., Introduction to Markov chain Monte Carlo, In
-- Handbook of Markov Chain Monte Carlo (pp. 45) (2011). CRC press.
module Mcmc.Algorithm.MHG
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
import Data.Time
import Mcmc.Algorithm
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Save
import Mcmc.Chain.Trace
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log
import System.Random.MWC
import Text.Printf
import Prelude hiding (cycle)

-- | The MHG algorithm.
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

-- NOTE: IO is required because the trace is mutable.

-- | Initialize an MHG algorithm.
mhg ::
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  TraceLength ->
  InitialState a ->
  GenIO ->
  IO (MHG a)
mhg pr lh cc mn trLen i0 g = do
  -- The trace is a mutable vector and the mutable state needs to be handled by
  -- a monad.
  tr <- replicateT traceLength l0
  return $ MHG $ Chain 0 l0 0 tr ac g 0 pr lh cc mn
  where
    l0 = Link i0 (pr i0) (lh i0)
    ac = emptyA $ ccProposals cc
    batchMonitorSizes = map getMonitorBatchSize $ mBatches mn
    minimumTraceLength = case trLen of
      TraceAuto -> 1
      TraceMinimum n -> n
    traceLength = maximum $ minimumTraceLength : batchMonitorSizes

mhgFn :: AnalysisName -> FilePath
mhgFn (AnalysisName nm) = nm ++ ".mcmc.mhg"

-- | Save an MHG algorithm.
mhgSave ::
  ToJSON a =>
  AnalysisName ->
  MHG a ->
  IO ()
mhgSave nm (MHG c) = do
  savedChain <- toSavedChain c
  BL.writeFile (mhgFn nm) $ compress $ encode savedChain

-- | Load an MHG algorithm.
--
-- See 'Mcmc.Mcmc.mcmcContinue'.
mhgLoad ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  AnalysisName ->
  IO (MHG a)
mhgLoad pr lh cc mn nm = do
  savedChain <- eitherDecode . decompress <$> BL.readFile (mhgFn nm)
  chain <- either error (fromSavedChain pr lh cc mn) savedChain
  return $ MHG chain

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

mhgPush :: MHG a -> IO (MHG a)
mhgPush (MHG c) = do
  t' <- pushT i t
  return $ MHG c {trace = t', iteration = succ n}
  where
    i = link c
    t = trace c
    n = iteration c

-- Ignore the number of capabilities. I have tried a lot of stuff, but the MHG
-- algorithm is just inherently sequential. Parallelization can be achieved by
-- having parallel prior and/or likelihood functions, or by using algorithms
-- running parallel chains such as 'MC3'.
mhgIterate :: ParallelizationMode -> MHG a -> IO (MHG a)
mhgIterate _ a = do
  ps <- orderProposals cc g
  a' <- foldM mhgPropose a ps
  mhgPush a'
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

mhgOpenMonitors :: AnalysisName -> ExecutionMode -> MHG a -> IO (MHG a)
mhgOpenMonitors nm em (MHG c) = do
  m' <- mOpen pre suf em m
  return $ MHG c {monitor = m'}
  where
    m = monitor c
    pre = fromAnalysisName nm
    suf = printf "%02d" $ chainId c

mhgExecuteMonitors ::
  Verbosity ->
  -- Starting time.
  UTCTime ->
  -- Total number of iterations.
  Int ->
  MHG a ->
  IO (Maybe BL.ByteString)
mhgExecuteMonitors vb t0 iTotal (MHG c) = mExec vb i i0 t0 tr iTotal m
  where
    i = iteration c
    i0 = start c
    tr = trace c
    m = monitor c

mhgStdMonitorHeader :: MHG a -> BL.ByteString
mhgStdMonitorHeader (MHG c) = msHeader (mStdOut $ monitor c)

mhgCloseMonitors :: MHG a -> IO (MHG a)
mhgCloseMonitors (MHG c) = do
  m' <- mClose m
  return $ MHG $ c {monitor = m'}
  where
    m = monitor c
