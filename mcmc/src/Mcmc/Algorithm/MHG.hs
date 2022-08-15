{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Algorithm.MHG
-- Description :  Metropolis-Hastings-Green algorithm
-- Copyright   :  2021 Dominik Schrempf
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
import Data.Maybe
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
import System.Random.Stateful
import Text.Printf
import Prelude hiding (cycle)

-- | The MHG algorithm.
newtype MHG a = MHG {fromMHG :: Chain a}

instance ToJSON a => Algorithm (MHG a) where
  aName = const "Metropolis-Hastings-Green (MHG)"
  aIteration = iteration . fromMHG
  aIsInvalidState = mhgIsInvalidState
  aIterate = mhgIterate
  aAutoTune = mhgAutoTune
  aResetAcceptance = mhgResetAcceptance
  aCleanAfterBurnIn = mhgCleanAfterBurnIn
  aSummarizeCycle = mhgSummarizeCycle
  aOpenMonitors = mhgOpenMonitors
  aExecuteMonitors = mhgExecuteMonitors
  aStdMonitorHeader = mhgStdMonitorHeader
  aCloseMonitors = mhgCloseMonitors
  aSave = mhgSave

-- Calculate required length of trace. The length may be larger during burn in,
-- because the tuners of some proposals (e.g., HMC, NUTS) require the states of
-- the last tuning interval.
getTraceLength ::
  Maybe BurnInSettings ->
  TraceLength ->
  Monitor a ->
  Cycle a ->
  Int
getTraceLength burnIn tl mn cc = maximum $ minimumTraceLength : bi : batchMonitorSizes
  where
    batchMonitorSizes = map getMonitorBatchSize $ mBatches mn
    minimumTraceLength = case tl of
      TraceAuto -> 1
      TraceMinimum n -> n
    bi = case (ccRequireTrace cc, burnIn) of
      (True, Just (BurnInWithAutoTuning _ n)) -> n
      (True, Just (BurnInWithCustomAutoTuning ns ms)) -> max (maximum $ 0 : ns) (maximum $ 0 : ms)
      _ -> 0

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
  StdGen ->
  IO (MHG a)
mhg s pr lh cc mn i0 g = do
  -- The trace is a mutable vector and the mutable state needs to be handled by
  -- a monad.
  tr <- replicateT tl l0
  gm <- newIOGenM g
  return $ MHG $ Chain Nothing l0 0 tr ac gm 0 pr lh cc mn
  where
    l0 = Link i0 (pr i0) (lh i0)
    ac = emptyA $ ccProposals cc
    tl = getTraceLength (Just $ sBurnIn s) (sTraceLength s) mn cc

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
-- Also create a backup of the save.
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

-- | Like 'mhgLoad' but do not perform sanity checks.
--
-- Also create a backup of the save.
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
  -- copyFile fn fnBak
  savedChain <- eitherDecode . decompress <$> BL.readFile fn
  chain <- either error (f pr lh cc mn) savedChain
  return $ MHG chain
  where
    -- fnBak = mhgFn $ AnalysisName $ (fromAnalysisName nm ++ ".bak")
    fn = mhgFn nm

-- | MHG ratios are stored in log domain.
type MHGRatio = Log Double

-- The MHG ratio. This implementation has the following properties:
--
-- - The kernel ratio and the Jacobian are checked carefully and should be
-- - strictly positive, finite numbers.
--
-- - The ratio is 'Infinity' if fX is zero. In this case, the proposal is always
--   accepted.
--
-- - The ratio is 'NaN' if fY and fX are zero. In this case, the proposal is
--   always rejected.
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
-- [1] Handbook of Markov chain Monte Carlo (2011), CRC press.
--
-- [2] Dellaportas, P., & Roberts, G. O., An introduction to MCMC, Lecture Notes
-- in Statistics, (), 1–41 (2003).
-- http://dx.doi.org/10.1007/978-0-387-21811-3_1.
mhgRatio :: Posterior -> Posterior -> KernelRatio -> Jacobian -> MHGRatio
mhgRatio fX fY q j
  | q == 0.0 = error "mhgRatio: Kernel ratio is negative infinity. Use 'ForceReject'."
  | q == 1.0 / 0.0 = error "mhgRatio: Kernel ratio is infinity. Use 'ForceAccept'."
  | q == 0.0 / 0.0 = error "mhgRatio: Kernel ratio is NaN."
  | j == 0.0 = error "mhgRatio: Jacobian is negative infinity. Use 'ForceReject'."
  | j == 1.0 / 0.0 = error "mhgRatio: Jacobian is infinity. Use 'ForceAccept'."
  | j == 0.0 / 0.0 = error "mhgRatio: Jacobian is NaN."
  | otherwise = fY / fX * q * j
{-# INLINE mhgRatio #-}

-- | Accept or reject a proposal with given MHG ratio?
mhgAccept :: MHGRatio -> IOGenM StdGen -> IO Bool
mhgAccept r g
  | ln r >= 0.0 = return True
  | otherwise = do
      b <- uniformRM (0, 1) g
      return $ b < exp (ln r)

mhgPropose :: MHG a -> Proposal a -> IO (MHG a)
mhgPropose (MHG c) p = do
  -- 1. Sample new state.
  !(pres, mcs) <- liftIO $ s x g
  -- 2. Define new prior and likelihood calculation functions. Avoid actual
  -- calculation of the values.
  --
  -- Most often, parallelization is not helpful, because the prior and
  -- likelihood functions are too fast; see
  -- https://stackoverflow.com/a/46603680/3536806.
  let calcPrLh y = (pF y, lF y) `using` parTuple2 rdeepseq rdeepseq
      accept y pr lh =
        let !ac' = case mcs of
              Nothing -> pushAccept p ac
              Just cs -> pushAcceptanceCounts p cs ac
         in pure $ MHG $ c {link = Link y pr lh, acceptance = ac'}
      reject =
        let !ac' = case mcs of
              Nothing -> pushReject p ac
              Just cs -> pushAcceptanceCounts p cs ac
         in pure $ MHG $ c {acceptance = ac'}
  -- 3. Accept or reject.
  --
  -- 3a. When rejection is inevitable, avoid calculation of the prior, the
  -- likelihood and the MHG ratio.
  case pres of
    ForceReject -> reject
    ForceAccept y -> let (pY, lY) = calcPrLh y in accept y pY lY
    (Propose y q j) ->
      if q <= 0.0 || j <= 0.0
        then reject
        else do
          -- 3b. Calculate Metropolis-Hastings-Green ratio.
          let (pY, lY) = calcPrLh y
              !r = mhgRatio (pX * lX) (pY * lY) q j
          isAccept <- mhgAccept r g
          if isAccept
            then accept y pY lY
            else reject
  where
    s = prFunction p
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
mhgIsInvalidState :: MHG a -> Bool
mhgIsInvalidState a = checkSoft p || check l || check (p * l)
  where
    x = link $ fromMHG a
    p = prior x
    l = likelihood x
    check v = let v' = ln v in isNaN v' || isInfinite v' || v' == 0
    checkSoft v = let v' = ln v in isNaN v' || isInfinite v'

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

mhgAutoTune :: TuningType -> Int -> MHG a -> IO (MHG a)
mhgAutoTune b n (MHG c) = do
  mxs <-
    if ccRequireTrace cc
      then Just . VB.map state <$> takeT n tr
      else pure Nothing
  return $ MHG $ c {cycle = autoTuneCycle b ac mxs cc}
  where
    ac = acceptance c
    cc = cycle c
    tr = trace c

mhgResetAcceptance :: MHG a -> MHG a
mhgResetAcceptance (MHG c) = MHG $ c {acceptance = resetA ac}
  where
    ac = acceptance c

mhgCleanAfterBurnIn :: TraceLength -> MHG a -> IO (MHG a)
mhgCleanAfterBurnIn tl (MHG c) = do
  xs <- takeT l tr
  tr' <- fromVectorT xs
  let c' = c {trace = tr'}
  pure $ MHG c'
  where
    mn = monitor c
    cc = cycle c
    tr = trace c
    l = getTraceLength Nothing tl mn cc

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
