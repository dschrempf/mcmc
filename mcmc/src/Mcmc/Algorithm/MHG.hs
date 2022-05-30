{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Algorithm.MHG
-- Description :  Metropolis-Hastings-Green algorithm
-- Copyright   :  (c) Dominik Schrempf 2021
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
    mhgLoadUnsafe,
    MHGRatio,
    mhgAccept,
  )
where

import Codec.Compression.GZip
import Control.Monad
import Control.Monad.IO.Class
import Control.Parallel.Strategies
import Data.Aeson
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Time
import qualified Data.Vector as VB
import Mcmc.Acceptance
import Mcmc.Algorithm
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Save
import Mcmc.Chain.Trace
import Mcmc.Cycle
import Mcmc.Likelihood
import Mcmc.Monitor
import Mcmc.Posterior
import Mcmc.Prior hiding (uniform)
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
  aIsInValidState = mhgIsInValidState
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
--
-- NOTE: Computation in the 'IO' Monad is necessary because the trace is
-- mutable.
mhg ::
  Settings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  InitialState a ->
  GenIO ->
  IO (MHG a)
mhg s pr lh cc mn i0 g = do
  -- The trace is a mutable vector and the mutable state needs to be handled by
  -- a monad.
  tr <- replicateT traceLength l0
  return $ MHG $ Chain Nothing l0 0 tr ac g 0 pr lh cc mn
  where
    l0 = Link i0 (pr i0) (lh i0)
    ac = emptyA $ ccProposals cc
    batchMonitorSizes = map getMonitorBatchSize $ mBatches mn
    minimumTraceLength = case sTraceLength s of
      TraceAuto -> 1
      TraceMinimum n -> n
    bi = case sBurnIn s of
      BurnInWithAutoTuning _ n -> n
      BurnInWithCustomAutoTuning ns ms -> max (maximum $ 0 : ns) (maximum $ 0 : ms)
      _ -> 0
    traceLength = maximum $ minimumTraceLength : bi : batchMonitorSizes

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
mhgLoad = mhgLoadWith fromSavedChain

-- | See 'mhgLoad' but do not perform sanity checks.
--
-- Useful when restarting a run with changed prior function, likelihood function
-- or proposals. Use with care!
mhgLoadUnsafe ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  AnalysisName ->
  IO (MHG a)
mhgLoadUnsafe = mhgLoadWith fromSavedChainUnsafe

-- Nice type :-).
mhgLoadWith ::
  FromJSON a =>
  (PriorFunction a -> LikelihoodFunction a -> Cycle a -> Monitor a -> SavedChain a -> IO (Chain a)) ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  AnalysisName ->
  IO (MHG a)
mhgLoadWith f pr lh cc mn nm = do
  savedChain <- eitherDecode . decompress <$> BL.readFile (mhgFn nm)
  chain <- either error (f pr lh cc mn) savedChain
  return $ MHG chain

-- | MHG ratios are stored in log domain.
type MHGRatio = Log Double

-- The MHG ratio. This implementation has the following properties:
--
-- - The ratio is 'Infinity' if fX is zero. In this case, the proposal is always
--   accepted.
--
-- - The ratio 'NaN' if (fY or q or j) and fX are zero. In this case, the
--   proposal is always rejected.
--
-- This means that a chain in a state with posterior probability zero (fX=0) can
-- only move if a state with non-zero posterior probability is proposed.
-- Otherwise it is stuck. Therefore, I print a warning when the posterior
-- probability is zero in the beginning of the MCMC run. This is probably not
-- the best behavior, but see below.
--
-- There is a discrepancy between authors saying that one should (a) always
-- accept the new state when the current posterior is zero (Chapter 4 of [1],
-- [2]), or (b) almost surely reject the proposal when either fY or q are zero
-- (Chapter 1 of [1]).
--
-- Since I trust the author of Chapter 1 (Charles Geyer) I choose to follow
-- option (b). However, Option (a) is more user-friendly.
--
-- [1] Handbook of markov chain monte carlo (2011), CRC press.
--
-- [2] Dellaportas, P., & Roberts, G. O., An introduction to mcmc, Lecture Notes
-- in Statistics, (), 1–41 (2003).
-- http://dx.doi.org/10.1007/978-0-387-21811-3_1.
mhgRatio :: Posterior -> Posterior -> KernelRatio -> Jacobian -> MHGRatio
-- q = qYX / qXY * jXY; see 'ProposalSimple'.
-- j = Jacobian.
mhgRatio fX fY q j = fY / fX * q * j
{-# INLINE mhgRatio #-}

-- | Accept or reject a proposal with given MHG ratio?
mhgAccept :: MHGRatio -> GenIO -> IO Bool
mhgAccept r g
  | ln r >= 0.0 = return True
  | otherwise = do
      b <- uniform g
      return $ b < exp (ln r)

mhgPropose :: MHG a -> Proposal a -> IO (MHG a)
mhgPropose (MHG c) p = do
  -- 1. Sample new state.
  (!y, !q, !j, mpy, mly) <- liftIO $ s x g
  -- 2. Calculate Metropolis-Hastings-Green ratio.
  let (pY, lY) = case (mpy, mly) of
        -- Most often, parallelization is not helpful, because the prior and
        -- likelihood functions are too fast; see
        -- https://stackoverflow.com/a/46603680/3536806.
        (Nothing, Nothing) -> (pF y, lF y) `using` parTuple2 rdeepseq rdeepseq
        (Just py, Nothing) -> (py, lF y)
        (Nothing, Just ly) -> (pF y, ly)
        (Just py, Just ly) -> (py, ly)
  -- let !pY = pF y
  --     !lY = lF y
  let !r = mhgRatio (pX * lX) (pY * lY) q j
  -- 3. Accept or reject.
  accept <- mhgAccept r g
  if accept
    then do
      let !ac' = pushA p True ac
      return $ MHG $ c {link = Link y pY lY, acceptance = ac'}
    else do
      let !ac' = pushA p False ac
      return $ MHG $ c {acceptance = pushA p False ac'}
  where
    s = prSimple p
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

-- Check if the current state is invalid.
--
-- At the moment this just checks whether the prior, likelihood, or posterior
-- are NaN or infinite.
mhgIsInValidState :: MHG a -> Bool
mhgIsInValidState a = check p || check l || check (p * l)
  where
    x = link $ fromMHG a
    p = prior x
    l = likelihood x
    check v = let v' = ln v in isNaN v' || isInfinite v' || v' == 0

-- Ignore the number of capabilities. I have tried a lot of stuff, but the MHG
-- algorithm is just inherently sequential. Parallelization can be achieved by
-- having parallel prior and/or likelihood functions, or by using algorithms
-- running parallel chains such as 'MC3'.
mhgIterate :: IterationMode -> ParallelizationMode -> MHG a -> IO (MHG a)
mhgIterate m _ a = do
  ps <- prepareProposals m cc g
  a' <- foldM mhgPropose a ps
  mhgPush a'
  where
    c = fromMHG a
    cc = cycle c
    g = generator c

mhgAutoTune :: Int -> MHG a -> IO (MHG a)
mhgAutoTune n (MHG c) = do
  tr <- VB.map state <$> takeT n (trace c)
  return $ MHG $ c {cycle = autoTuneCycle ac tr cc}
  where
    ac = acceptance c
    cc = cycle c

mhgResetAcceptance :: MHG a -> MHG a
mhgResetAcceptance (MHG c) = MHG $ c {acceptance = resetA ac}
  where
    ac = acceptance c

mhgSummarizeCycle :: IterationMode -> MHG a -> BL.ByteString
mhgSummarizeCycle m (MHG c) = summarizeCycle m ac cc
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
    suf = maybe "" (printf "%02d") $ chainId c

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
