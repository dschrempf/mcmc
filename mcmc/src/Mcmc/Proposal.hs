{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- TODO: Proposals on simplices: SimplexElementScale (?).

-- |
-- Module      :  Mcmc.Proposal
-- Description :  Proposals and cycles
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 13:42:53 2020.
module Mcmc.Proposal
  ( -- * Proposal
    Proposal (..),
    (@~),
    ProposalSimple (..),
    Tuner (tParam, tFunc),
    createProposal,
    tune,

    -- * Cycle
    Order (..),
    Cycle (ccProposals),
    fromList,
    setOrder,
    getNIterations,
    tuneCycle,
    autotuneCycle,
    summarizeCycle,

    -- * Acceptance
    Acceptance (fromAcceptance),
    emptyA,
    pushA,
    resetA,
    transformKeysA,
    acceptanceRatios,
  )
where

import Data.Aeson
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Default
import qualified Data.Double.Conversion.ByteString as BC
import Data.Function
import Data.List
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as M
import Data.Maybe
import Lens.Micro
import Mcmc.Internal.ByteString
import Mcmc.Tools.Shuffle
import Numeric.Log hiding (sum)
import System.Random.MWC

-- | A 'Proposal' is an instruction about how the Markov chain will traverse the
-- state space @a@. Essentially, it is a probability mass or probability density
-- conditioned on the current state (i.e., a kernel).
--
-- A 'Proposal' may be tuneable in that it contains information about how to enlarge
-- or shrink the step size to tune the acceptance ratio.
data Proposal a = Proposal
  { -- | Name (no proposals with the same name are allowed in a 'Cycle').
    pName :: String,
    -- | The weight determines how often a 'Proposal' is executed per iteration of
    -- the Markov chain.
    pWeight :: Int,
    -- | Simple proposal without name, weight, and tuning information.
    pSimple :: ProposalSimple a,
    -- | Tuning is disabled if set to 'Nothing'.
    pTuner :: Maybe (Tuner a)
  }

instance Show (Proposal a) where
  show m = show $ pName m

instance Eq (Proposal a) where
  m == n = pName m == pName n

instance Ord (Proposal a) where
  compare = compare `on` pName

-- | Convert a proposal from one data type to another using a lens.
--
-- For example:
--
-- @
-- scaleFirstEntryOfTuple = _1 @~ scale
-- @
(@~) :: Lens' b a -> Proposal a -> Proposal b
(@~) l (Proposal n w s t) = Proposal n w (convertS l s) (convertT l <$> t)

-- | Simple proposal without tuning information.
--
-- Instruction about randomly moving from the current state to a new state,
-- given some source of randomness.
--
-- In order to calculate the Metropolis-Hastings ratio, we need to know the
-- ratio of the backward to forward kernels (i.e., the probability masses or
-- probability densities). For unbiased proposals, this ratio is 1.0.
newtype ProposalSimple a = ProposalSimple
  { pSample :: a -> GenIO -> IO (a, Log Double)
  }

convertS :: Lens' b a -> ProposalSimple a -> ProposalSimple b
convertS l (ProposalSimple s) = ProposalSimple s'
  where
    s' v g = do
      (x', r) <- s (v ^. l) g
      return (set l x' v, r)

-- | Tune the acceptance ratio of a 'Proposal'; see 'tune', or 'autotuneCycle'.
data Tuner a = Tuner
  { tParam :: Double,
    tFunc :: Double -> ProposalSimple a
  }

convertT :: Lens' b a -> Tuner a -> Tuner b
convertT l (Tuner p f) = Tuner p f'
  where
    f' x = convertS l $ f x

-- | Create a possibly tuneable proposal.
createProposal ::
  -- | Name.
  String ->
  -- | Weight.
  Int ->
  -- | Function creating a simple proposal for a given tuning parameter. The
  -- larger the tuning parameter, the larger the proposal (and the lower the
  -- expected acceptance ratio), and vice versa.
  (Double -> ProposalSimple a) ->
  -- | Activate tuning?
  Bool ->
  Proposal a
createProposal n w f True = Proposal n w (f 1.0) (Just $ Tuner 1.0 f)
createProposal n w f False = Proposal n w (f 1.0) Nothing

-- Minimal tuning parameter; subject to change.
tuningParamMin :: Double
tuningParamMin = 1e-12

-- | Tune a 'Proposal'. Return 'Nothing' if 'Proposal' is not tuneable. If the parameter
--   @dt@ is larger than 1.0, the 'Proposal' is enlarged, if @0<dt<1.0@, it is
--   shrunk. Negative tuning parameters are not allowed.
tune :: Double -> Proposal a -> Maybe (Proposal a)
tune dt m
  | dt <= 0 = error $ "tune: Tuning parameter not positive: " <> show dt <> "."
  | otherwise = do
    (Tuner t f) <- pTuner m
    -- Ensure that the tuning parameter is not too small.
    let t' = max tuningParamMin (t * dt)
    return $ m {pSimple = f t', pTuner = Just $ Tuner t' f}

-- The desired acceptance ratio 0.44 is optimal for one-dimensional proposals;
-- one could also store the affected number of dimensions with the proposal and
-- tune towards an acceptance ratio accounting for the number of dimensions.
--
-- The optimal ratios seem to be:
-- - One dimension: 0.44 (numerical result).
-- - Five and more dimensions: 0.234 seems to be a good value (numerical result).
-- - Infinite dimensions: 0.234 (theorem for specific target distributions).
-- See Handbook of Markov chain Monte Carlo, chapter 4.
ratioOpt :: Double
ratioOpt = 0.44

-- Warn if acceptance ratio is lower.
ratioMin :: Double
ratioMin = 0.1

-- Warn if acceptance ratio is larger.
ratioMax :: Double
ratioMax = 0.9

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

instance Default Order where def = RandomO

-- | In brief, a 'Cycle' is a list of proposals.
--
-- The state of the Markov chain will be logged only after all 'Proposal's in
-- the 'Cycle' have been completed, and the iteration counter will be increased
-- by one. The order in which the 'Proposal's are executed is specified by
-- 'Order'. The default is 'RandomO'.
--
-- __Proposals must have unique names__, so that they can be identified.
data Cycle a = Cycle
  { ccProposals :: [Proposal a],
    ccOrder :: Order
  }

-- | Create a 'Cycle' from a list of 'Proposal's.
fromList :: [Proposal a] -> Cycle a
fromList [] =
  error "fromList: Received an empty list but cannot create an empty Cycle."
fromList xs =
  if length (nub nms) == length nms
    then Cycle xs def
    else error "fromList: Proposals don't have unique names."
  where
    nms = map pName xs

-- | Set the order of 'Proposal's in a 'Cycle'.
setOrder :: Order -> Cycle a -> Cycle a
setOrder o c = c {ccOrder = o}

-- | Replicate 'Proposal's according to their weights and possibly shuffle them.
getNIterations :: Cycle a -> Int -> GenIO -> IO [[Proposal a]]
getNIterations (Cycle xs o) n g = case o of
  RandomO -> shuffleN ps n g
  SequentialO -> return $ replicate n ps
  RandomReversibleO -> do
    psRs <- shuffleN ps n g
    return [psR ++ reverse psR | psR <- psRs]
  SequentialReversibleO -> return $ replicate n $ ps ++ reverse ps
  where
    !ps = concat [replicate (pWeight m) m | m <- xs]

-- | Tune 'Proposal's in the 'Cycle'. See 'tune'.
tuneCycle :: Map (Proposal a) Double -> Cycle a -> Cycle a
tuneCycle m c =
  if sort (M.keys m) == sort ps
    then c {ccProposals = map tuneF ps}
    else error "tuneCycle: Map contains proposals that are not in the cycle."
  where
    ps = ccProposals c
    tuneF p = case m M.!? p of
      Nothing -> p
      Just x -> fromMaybe p (tune x p)

-- | Calculate acceptance ratios and auto tune the 'Proposal's in the 'Cycle'. For
-- now, a 'Proposal' is enlarged when the acceptance ratio is above 0.44, and
-- shrunk otherwise. Do not change 'Proposal's that are not tuneable.
autotuneCycle :: Acceptance (Proposal a) -> Cycle a -> Cycle a
autotuneCycle a = tuneCycle (M.map (\x -> exp $ x - ratioOpt) $ acceptanceRatios a)

renderRow ::
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString
renderRow name weight nAccept nReject acceptRatio tuneParam manualAdjustment = "   " <> nm <> wt <> na <> nr <> ra <> tp <> mt
  where
    nm = alignLeft 30 name
    wt = alignRight 8 weight
    na = alignRight 15 nAccept
    nr = alignRight 15 nReject
    ra = alignRight 15 acceptRatio
    tp = alignRight 20 tuneParam
    mt = alignRight 30 manualAdjustment

proposalHeader :: BL.ByteString
proposalHeader =
  renderRow "Proposal" "Weight" "Accepted" "Rejected" "Ratio" "Tuning parameter" "Consider manual adjustment"

summarizeProposal :: Proposal a -> Maybe (Int, Int, Double) -> BL.ByteString
summarizeProposal m r =
  renderRow
    (BL.pack name)
    weight
    nAccept
    nReject
    acceptRatio
    tuneParamStr
    manualAdjustmentStr
  where
    name = pName m
    weight = BB.toLazyByteString $ BB.intDec $ pWeight m
    nAccept = BB.toLazyByteString $ maybe "" (BB.intDec . (^. _1)) r
    nReject = BB.toLazyByteString $ maybe "" (BB.intDec . (^. _2)) r
    acceptRatio = BL.fromStrict $ maybe "" (BC.toFixed 3 . (^. _3)) r
    tuneParamStr = BL.fromStrict $ maybe "" (BC.toFixed 3) (tParam <$> pTuner m)
    check v
      | v < ratioMin = "ratio too low"
      | v > ratioMax = "ratio too high"
      | otherwise = ""
    manualAdjustmentStr = BL.fromStrict $ maybe "" (check . (^. _3)) r

hLine :: BL.ByteString -> BL.ByteString
hLine s = "   " <> BL.replicate (BL.length s - 3) '-'

-- | Summarize the 'Proposal's in the 'Cycle'. Also report acceptance ratios.
summarizeCycle :: Acceptance (Proposal a) -> Cycle a -> BL.ByteString
summarizeCycle a c =
  BL.intercalate "\n" $
    [ "Summary of proposal(s) in cycle. " <> mpi <> " proposal(s) per iteration.",
      proposalHeader,
      hLine proposalHeader
    ]
      ++ [summarizeProposal m (ar m) | m <- ps]
      ++ [hLine proposalHeader]
  where
    ps = ccProposals c
    mpi = BB.toLazyByteString $ BB.intDec $ sum $ map pWeight ps
    ar m = acceptanceRatio m a

-- | For each key @k@, store the number of accepted and rejected proposals.
newtype Acceptance k = Acceptance {fromAcceptance :: Map k (Int, Int)}

instance ToJSONKey k => ToJSON (Acceptance k) where
  toJSON (Acceptance m) = toJSON m
  toEncoding (Acceptance m) = toEncoding m

instance (Ord k, FromJSONKey k) => FromJSON (Acceptance k) where
  parseJSON v = Acceptance <$> parseJSON v

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
emptyA :: Ord k => [k] -> Acceptance k
emptyA ks = Acceptance $ M.fromList [(k, (0, 0)) | k <- ks]

-- | For key @k@, prepend an accepted (True) or rejected (False) proposal.
pushA :: (Ord k, Show k) => k -> Bool -> Acceptance k -> Acceptance k
pushA k True = Acceptance . M.adjust (\(a, r) -> (succ a, r)) k . fromAcceptance
pushA k False = Acceptance . M.adjust (\(a, r) -> (a, succ r)) k . fromAcceptance
{-# INLINEABLE pushA #-}

-- | Reset acceptance storage.
resetA :: Ord k => Acceptance k -> Acceptance k
resetA = emptyA . M.keys . fromAcceptance

transformKeys :: (Ord k1, Ord k2) => [k1] -> [k2] -> Map k1 v -> Map k2 v
transformKeys ks1 ks2 m = foldl' insrt M.empty $ zip ks1 ks2
  where
    insrt m' (k1, k2) = M.insert k2 (m M.! k1) m'

-- | Transform keys using the given lists. Keys not provided will not be present
-- in the new 'Acceptance' variable.
transformKeysA :: (Ord k1, Ord k2) => [k1] -> [k2] -> Acceptance k1 -> Acceptance k2
transformKeysA ks1 ks2 = Acceptance . transformKeys ks1 ks2 . fromAcceptance

-- | Acceptance counts and ratio for a specific proposal.
acceptanceRatio :: (Show k, Ord k) => k -> Acceptance k -> Maybe (Int, Int, Double)
acceptanceRatio k a = case fromAcceptance a M.!? k of
  Just (0, 0) -> Nothing
  Just (as, rs) -> Just (as, rs, fromIntegral as / fromIntegral (as + rs))
  Nothing -> error $ "acceptanceRatio: Key not found in map: " ++ show k ++ "."

-- | Acceptance ratios for all proposals.
acceptanceRatios :: Acceptance k -> Map k Double
acceptanceRatios = M.map (\(as, rs) -> fromIntegral as / fromIntegral (as + rs)) . fromAcceptance
