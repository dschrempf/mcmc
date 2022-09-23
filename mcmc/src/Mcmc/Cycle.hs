{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Cycle
-- Description :  A cycle is a list of proposals
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Thu Jul  8 17:56:03 2021.
module Mcmc.Cycle
  ( -- * Cycles
    Order (..),
    Cycle (ccProposals, ccRequireTrace, ccHasIntermediateTuners),
    cycleFromList,
    setOrder,
    IterationMode (..),
    prepareProposals,
    autoTuneCycle,

    -- * Output
    summarizeCycle,
  )
where

import Control.Applicative
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List
import qualified Data.Map.Strict as M
import Data.Maybe
import qualified Data.Vector as VB
import Mcmc.Acceptance
import Mcmc.Internal.Shuffle
import Mcmc.Proposal
import System.Random.Stateful

-- | Define the order in which 'Proposal's are executed in a 'Cycle'. The total
-- number of 'Proposal's per 'Cycle' may differ between 'Order's (e.g., compare
-- 'RandomO' and 'RandomReversibleO').
data Order
  = -- | Shuffle the 'Proposal's in the 'Cycle'. The 'Proposal's are replicated
    -- according to their weights and executed in random order. If a 'Proposal' has
    -- weight @w@, it is executed exactly @w@ times per iteration.
    RandomO
  | -- | The 'Proposal's are executed sequentially, in the order they appear in the
    -- 'Cycle'. 'Proposal's with weight @w>1@ are repeated immediately @w@ times
    -- (and not appended to the end of the list).
    SequentialO
  | -- | Similar to 'RandomO'. However, a reversed copy of the list of
    --  shuffled 'Proposal's is appended such that the resulting Markov chain is
    --  reversible.
    --  Note: the total number of 'Proposal's executed per cycle is twice the number
    --  of 'RandomO'.
    RandomReversibleO
  | -- | Similar to 'SequentialO'. However, a reversed copy of the list of
    -- sequentially ordered 'Proposal's is appended such that the resulting Markov
    -- chain is reversible.
    SequentialReversibleO
  deriving (Eq, Show)

-- Describe the order.
describeOrder :: Order -> BL.ByteString
describeOrder RandomO = "The proposals are executed in random order."
describeOrder SequentialO = "The proposals are executed sequentially."
describeOrder RandomReversibleO =
  BL.intercalate
    "\n"
    [ describeOrder RandomO,
      "A reversed copy of the shuffled proposals is appended to ensure reversibility."
    ]
describeOrder SequentialReversibleO =
  BL.intercalate
    "\n"
    [ describeOrder SequentialO,
      "A reversed copy of the sequential proposals is appended to ensure reversibility."
    ]

-- | In brief, a 'Cycle' is a list of proposals.
--
-- The state of the Markov chain will be logged only after all 'Proposal's in
-- the 'Cycle' have been completed, and the iteration counter will be increased
-- by one. The order in which the 'Proposal's are executed is specified by
-- 'Order'. The default is 'RandomO'.
--
-- No proposals with the same name and description are allowed in a 'Cycle', so
-- that they can be uniquely identified.
data Cycle a = Cycle
  { ccProposals :: [Proposal a],
    ccOrder :: Order,
    -- | Does the cycle require the trace when auto tuning? See 'tRequireTrace'.
    ccRequireTrace :: Bool,
    -- | Does the cycle include proposals that can be tuned every iterations?
    -- See 'tSuitableForIntermediateTuning'.
    ccHasIntermediateTuners :: Bool
  }

-- | Create a 'Cycle' from a list of 'Proposal's; use 'RandomO', but see 'setOrder'.
cycleFromList :: [Proposal a] -> Cycle a
cycleFromList [] =
  error "cycleFromList: Received an empty list but cannot create an empty Cycle."
cycleFromList xs =
  if length uniqueXs == length xs
    then Cycle xs RandomO (any needsTrace xs) (any isIntermediate xs)
    else error $ "\n" ++ msg ++ "cycleFromList: Proposals are not unique."
  where
    uniqueXs = nub xs
    removedXs = xs \\ uniqueXs
    removedNames = map (show . prName) removedXs
    removedDescriptions = map (show . prDescription) removedXs
    removedMsgs = zipWith (\n d -> n ++ " " ++ d) removedNames removedDescriptions
    msg = unlines removedMsgs
    needsTrace p = maybe False tRequireTrace (prTuner p)
    isIntermediate p = maybe False tSuitableForIntermediateTuning (prTuner p)

-- | Set the order of 'Proposal's in a 'Cycle'.
setOrder :: Order -> Cycle a -> Cycle a
setOrder o c = c {ccOrder = o}

-- | Use all proposals, or use fast proposals only?
data IterationMode = AllProposals | FastProposals
  deriving (Eq)

-- | Replicate 'Proposal's according to their weights and possibly shuffle them.
prepareProposals :: StatefulGen g m => IterationMode -> Cycle a -> g -> m [Proposal a]
prepareProposals m (Cycle xs o _ _) g =
  if null ps
    then
      let msg = case m of
            FastProposals -> "no fast proposals found"
            AllProposals -> "no proposals found"
       in error $ "prepareProposals: " <> msg
    else case o of
      RandomO -> shuffle ps g
      SequentialO -> return ps
      RandomReversibleO -> do
        psR <- shuffle ps g
        return $ psR ++ reverse psR
      SequentialReversibleO -> return $ ps ++ reverse ps
  where
    !ps =
      concat
        [ replicate (fromPWeight $ prWeight p) p
          | p <- xs,
            case m of
              AllProposals -> True
              -- Only use proposal if it is fast.
              FastProposals -> prSpeed p == PFast
        ]

-- The number of proposals depends on the order.
getNProposalsPerCycle :: IterationMode -> Cycle a -> Int
getNProposalsPerCycle m (Cycle xs o _ _) = case o of
  RandomO -> once
  SequentialO -> once
  RandomReversibleO -> 2 * once
  SequentialReversibleO -> 2 * once
  where
    xs' = case m of
      AllProposals -> xs
      FastProposals -> filter (\x -> prSpeed x == PFast) xs
    once = sum $ map (fromPWeight . prWeight) xs'

-- See 'tuneWithTuningParameters' and 'Tuner'.
tuneWithChainParameters ::
  TuningType ->
  Maybe AcceptanceRate ->
  Maybe (VB.Vector a) ->
  Proposal a ->
  Either String (Proposal a)
tuneWithChainParameters tt mar mxs p = case prTuner p of
  Nothing -> Right p
  Just (Tuner t ts rt it fT _) -> case (tt, it, prSpeed p) of
    (IntermediateTuningFastProposalsOnly, True, PFast) -> tuneIntermediate
    (IntermediateTuningAllProposals, True, _) -> tuneIntermediate
    (NormalTuningFastProposalsOnly, _, PFast) -> tuneNormally
    (NormalTuningAllProposals, _, _) -> tuneNormally
    (LastTuningFastProposalsOnly, _, _) -> tuneNormally
    (LastTuningAllProposals, _, _) -> tuneNormally
    _ -> Right p
    where
      hasTrace = isJust mxs
      err m = Left $ "tuneWithChainParameters: " <> m
      tuneIntermediate =
        if hasTrace
          then err "intermediate tuning but trace provided"
          else tune
      tuneNormally =
        if rt && not hasTrace
          then err "trace required"
          else tune
      tune =
        let (t', ts') = fT tt (prDimension p) mar mxs (t, ts)
         in tuneWithTuningParameters t' ts' p

-- (_, False, Just _) ->
-- (True, _, Just _) ->
-- (False, True, Nothing) ->
-- _ ->

-- | Calculate acceptance rates and auto tunes the 'Proposal's in the 'Cycle'.
--
-- Do not change 'Proposal's that are not tuneable.
autoTuneCycle :: TuningType -> Acceptances (Proposal a) -> Maybe (VB.Vector a) -> Cycle a -> Cycle a
autoTuneCycle tt a mxs c
  | isJust mxs && not (ccRequireTrace c) = err "trace provided but not required"
  | otherwise =
      if sort (M.keys $ fromAcceptances a) == sort ps
        then c {ccProposals = map tuneF ps}
        else err "proposals in map and cycle do not match"
  where
    err msg = error $ "autoTuneCycle: " <> msg
    ps = ccProposals c
    tuneF p =
      let (_, _, mar, mtr) = acceptanceRate p a
          -- Favor the expected rate, if available.
          mr = mtr <|> mar
       in either error id $ tuneWithChainParameters tt mr mxs p

-- | Summarize the 'Proposal's in the 'Cycle'. Also report acceptance rates.
summarizeCycle :: IterationMode -> Acceptances (Proposal a) -> Cycle a -> BL.ByteString
summarizeCycle m a c =
  BL.intercalate "\n" $
    [ "Summary of proposal(s) in cycle.",
      nProposalsFullStr,
      describeOrder (ccOrder c),
      proposalHeader,
      proposalHLine
    ]
      ++ [ summarizeProposal
             (prName p)
             (prDescription p)
             (prWeight p)
             (tTuningParameter <$> prTuner p)
             (prDimension p)
             (ar p)
           | p <- ps
         ]
      ++ [proposalHLine]
  where
    ps = ccProposals c
    nProposals = getNProposalsPerCycle m c
    nProposalsStr = BB.toLazyByteString $ BB.intDec nProposals
    nProposalsFullStr = case nProposals of
      1 -> nProposalsStr <> " proposal is performed per iteration."
      _ -> nProposalsStr <> " proposals are performed per iterations."
    ar pr = acceptanceRate pr a
    proposalHLine = BL.replicate (BL.length proposalHeader) '-'
