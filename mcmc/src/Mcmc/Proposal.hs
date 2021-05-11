{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DerivingVia #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal
-- Description :  Proposals are instruction to move around the state space
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 13:42:53 2020.
module Mcmc.Proposal
  ( -- * Proposals and types
    PName (..),
    PDescription (..),
    PWeight (..),
    PDimension (..),
    Proposal (..),
    (@~),
    ProposalSimple,
    Tuner (..),
    Tune (..),
    createProposal,
    TuningParameter,
    tuningParameterMin,
    tuningParameterMax,
    tune,
    getOptimalRate,

    -- * Cycles
    Order (..),
    Cycle (ccProposals),
    cycleFromList,
    setOrder,
    prepareProposals,
    tuneCycle,
    autoTuneCycle,

    -- * Acceptance rates
    Acceptance (fromAcceptance),
    emptyA,
    pushA,
    resetA,
    transformKeysA,
    acceptanceRate,
    acceptanceRates,

    -- * Output
    proposalHeader,
    proposalHLine,
    summarizeProposal,
    summarizeCycle,
  )
where

import Control.DeepSeq
import Data.Aeson
import Data.Bifunctor
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Default
import qualified Data.Double.Conversion.ByteString as BC
import Data.Function
import Data.List
import qualified Data.Map.Strict as M
import Data.Maybe
import Lens.Micro
import Mcmc.Internal.ByteString
import Mcmc.Internal.Shuffle
import Numeric.Log hiding (sum)
import System.Random.MWC

-- | Proposal name.
newtype PName = PName {fromPName :: String}
  deriving (Show, Eq, Ord)
  deriving (Monoid, Semigroup) via String

-- | Proposal description.
newtype PDescription = PDescription {fromPDescription :: String}
  deriving (Show, Eq, Ord)

-- | The weight determines how often a 'Proposal' is executed per iteration of
-- the Markov chain.
newtype PWeight = PWeight {fromPWeight :: Int}
  deriving (Show, Eq, Ord)

-- | Proposal dimension.
--
-- The number of affected, independent parameters.
--
-- The optimal acceptance rate of low dimensional proposals is higher than for
-- high dimensional ones.
--
-- Optimal acceptance rates are still subject to controversies. As far as I
-- know, research has focused on random walk proposal with a multivariate normal
-- distribution of dimension @d@. In this case, the following acceptance rates
-- are desired:
--
-- - one dimension: 0.44 (numerical results);
--
-- - five and more dimensions: 0.234 (numerical results);
--
-- - infinite dimensions: 0.234 (theorem for specific target distributions).
--
-- See Handbook of Markov chain Monte Carlo, chapter 4.
--
-- Of course, many proposals may not be classical random walk proposals. For
-- example, the beta proposal on a simplex ('Mcmc.Proposal.Simplex.beta')
-- samples one new variable of the simplex from a beta distribution while
-- rescaling all other variables. What is the dimension of this proposal? I
-- don't know, but I set the dimension to 2. The reason is that if the dimension
-- of the simplex is 2, two variables are changed. If the dimension of the
-- simplex is high, one variable is changed substantially, while all others are
-- changed marginally.
--
-- Further, if a proposal changes a number of variables in the same way (and not
-- independently like in a random walk proposal), I still set the dimension of
-- the proposal to the number of variables changed.
--
-- Finally, I assume that proposals of unknown dimension have high dimension,
-- and use the optimal acceptance rate 0.234.
data PDimension = PDimension Int | PDimensionUnknown

-- | A 'Proposal' is an instruction about how the Markov chain will traverse the
-- state space @a@. Essentially, it is a probability mass or probability density
-- conditioned on the current state (i.e., a Markov kernel).
--
-- A 'Proposal' may be tuneable in that it contains information about how to enlarge
-- or shrink the step size to tune the acceptance rate.
data Proposal a = Proposal
  { -- | Name of the affected variable.
    pName :: PName,
    -- | Description of the proposal type and parameters.
    pDescription :: PDescription,
    -- | Dimension of the proposal. The dimension is used to calculate the
    -- optimal acceptance rate, and does not have to be exact.
    pDimension :: PDimension,
    -- | The weight determines how often a 'Proposal' is executed per iteration of
    -- the Markov chain.
    pWeight :: PWeight,
    -- | Simple proposal without name, weight, and tuning information.
    pSimple :: ProposalSimple a,
    -- | Tuning is disabled if set to 'Nothing'.
    pTuner :: Maybe (Tuner a)
  }

instance Eq (Proposal a) where
  m == n = pName m == pName n && pDescription m == pDescription n

instance Ord (Proposal a) where
  compare = compare `on` (\p -> (pDescription p, pName p, pWeight p))

-- | Convert a proposal from one data type to another using a lens.
--
-- For example:
--
-- @
-- scaleFirstEntryOfTuple = _1 @~ scale
-- @
(@~) :: Lens' b a -> Proposal a -> Proposal b
(@~) l (Proposal n r d w s t) = Proposal n r d w (convertProposalSimple l s) (convertTuner l <$> t)

-- | Simple proposal without tuning information.
--
-- Instruction about randomly moving from the current state to a new state,
-- given some source of randomness.
--
-- In order to calculate the Metropolis-Hastings-Green ratio, we need to know
-- the ratio of the backward to forward kernels (i.e., the probability masses or
-- probability densities) and the absolute value of the determinant of the
-- Jacobian matrix.
--
-- For unbiased proposals, these values are 1.0 such that
--
-- @
-- proposalSimpleUnbiased x g = return (x', 1.0, 1.0)
-- @
--
-- For biased proposals, the kernel ratio is qYX / qXY, where qXY is the
-- probability density to move from X to Y, and the absolute value of the
-- determinant of the Jacobian matrix differs from 1.0.
type ProposalSimple a = a -> GenIO -> IO (a, Log Double, Log Double)

convertProposalSimple :: Lens' b a -> ProposalSimple a -> ProposalSimple b
convertProposalSimple l s = s'
  where
    s' v g = do
      (x', r, j) <- s (v ^. l) g
      return (set l x' v, r, j)

-- | Tune the acceptance rate of a 'Proposal'; see 'tune', or 'autoTuneCycle'.
data Tuner a = Tuner
  { tParam :: TuningParameter,
    tFunc :: TuningParameter -> ProposalSimple a
  }

convertTuner :: Lens' b a -> Tuner a -> Tuner b
convertTuner l (Tuner p f) = Tuner p f'
  where
    f' x = convertProposalSimple l $ f x

-- | Tune the proposal?
data Tune = Tune | NoTune
  deriving (Show, Eq)

-- | Tuning parameter.
type TuningParameter = Double

-- | Create a tuneable proposal.
createProposal ::
  -- | Description of the proposal type and parameters.
  PDescription ->
  -- | Function creating a simple proposal for a given tuning parameter. The
  -- larger the tuning parameter, the larger the proposal (and the lower the
  -- expected acceptance rate), and vice versa.
  (TuningParameter -> ProposalSimple a) ->
  -- | Dimension.
  PDimension ->
  -- | Name.
  PName ->
  -- | Weight.
  PWeight ->
  -- | Activate tuning?
  Tune ->
  Proposal a
createProposal r f d n w Tune = Proposal n r d w (f 1.0) (Just $ Tuner 1.0 f)
createProposal r f d n w NoTune = Proposal n r d w (f 1.0) Nothing

-- | Minimal tuning parameter; @1e-12@, subject to change.
--
-- >>> tuningParameterMin
-- 1e-5
tuningParameterMin :: TuningParameter
tuningParameterMin = 1e-5

-- | Maximal tuning parameter; @1e12@, subject to change.
-- >>> tuningParameterMax
-- 1e3
tuningParameterMax :: TuningParameter
tuningParameterMax = 1e3

-- | Tune a 'Proposal'.
--
-- The size of the proposal is proportional to the tuning parameter which has a
-- positive lower bound of 'tuningParameterMin'.
--
-- The tuning function maps the current tuning parameter to a new one.
--
-- Return 'Nothing' if 'Proposal' is not tuneable.
tune :: (TuningParameter -> TuningParameter) -> Proposal a -> Maybe (Proposal a)
tune f m = do
  (Tuner t g) <- pTuner m
  -- Ensure that the tuning parameter is strictly positive and well bounded.
  let t' = max tuningParameterMin (f t)
      t'' = min tuningParameterMax t'
  return $ m {pSimple = g t'', pTuner = Just $ Tuner t'' g}

-- | See 'PDimension'.
getOptimalRate :: PDimension -> Double
getOptimalRate (PDimension n)
  | n <= 0 = error "getOptimalRate: Proposal dimension is zero or negative."
  | n == 1 = 0.44
  -- Use a linear interpolation with delta 0.0515.
  | n == 2 = 0.3885
  | n == 3 = 0.337
  | n == 4 = 0.2855
  | n >= 5 = 0.234
  | otherwise = error "getOptimalRate: Proposal dimension is not an integer?"
getOptimalRate PDimensionUnknown = 0.234

-- Warn if acceptance rate is lower.
rateMin :: Double
rateMin = 0.1

-- Warn if acceptance rate is larger.
rateMax :: Double
rateMax = 0.9

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
    ccOrder :: Order
  }

-- | Create a 'Cycle' from a list of 'Proposal's.
cycleFromList :: [Proposal a] -> Cycle a
cycleFromList [] =
  error "cycleFromList: Received an empty list but cannot create an empty Cycle."
cycleFromList xs =
  if length (nub xs) == length xs
    then Cycle xs def
    else error "cycleFromList: Proposals are not unique."

-- | Set the order of 'Proposal's in a 'Cycle'.
setOrder :: Order -> Cycle a -> Cycle a
setOrder o c = c {ccOrder = o}

-- | Replicate 'Proposal's according to their weights and possibly shuffle them.
prepareProposals :: Cycle a -> GenIO -> IO [Proposal a]
prepareProposals (Cycle xs o) g = case o of
  RandomO -> shuffle ps g
  SequentialO -> return ps
  RandomReversibleO -> do
    psR <- shuffle ps g
    return $ psR ++ reverse psR
  SequentialReversibleO -> return $ ps ++ reverse ps
  where
    !ps = concat [replicate (fromPWeight $ pWeight p) p | p <- xs]

-- The number of proposals depends on the order.
getNProposalsPerCycle :: Cycle a -> Int
getNProposalsPerCycle (Cycle xs o) = case o of
  RandomO -> once
  SequentialO -> once
  RandomReversibleO -> 2 * once
  SequentialReversibleO -> 2 * once
  where
    once = sum $ map (fromPWeight . pWeight) xs

-- | Tune 'Proposal's in the 'Cycle'. See 'tune'.
tuneCycle :: M.Map (Proposal a) (TuningParameter -> TuningParameter) -> Cycle a -> Cycle a
tuneCycle m c =
  if sort (M.keys m) == sort ps
    then c {ccProposals = map tuneF ps}
    else error "tuneCycle: Propoals in map and cycle do not match."
  where
    ps = ccProposals c
    tuneF p = case m M.!? p of
      Nothing -> p
      Just f -> fromMaybe p (tune f p)

-- | Calculate acceptance rates and auto tune the 'Proposal's in the 'Cycle'. For
-- now, a 'Proposal' is enlarged when the acceptance rate is above 0.44, and
-- shrunk otherwise. Do not change 'Proposal's that are not tuneable.
autoTuneCycle :: Acceptance (Proposal a) -> Cycle a -> Cycle a
autoTuneCycle a = tuneCycle (M.mapWithKey tuningF $ acceptanceRates a)
  where
    tuningF proposal mCurrentRate currentTuningParam = case mCurrentRate of
      Nothing -> currentTuningParam
      Just currentRate ->
        let optimalRate = getOptimalRate (pDimension proposal)
         in exp (2 * (currentRate - optimalRate)) * currentTuningParam

renderRow ::
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString ->
  BL.ByteString
renderRow name ptype weight nAccept nReject acceptRate optimalRate tuneParam manualAdjustment = nm <> pt <> wt <> na <> nr <> ra <> ro <> tp <> mt
  where
    nm = alignLeft 30 name
    pt = alignLeft 50 ptype
    wt = alignRight 8 weight
    na = alignRight 14 nAccept
    nr = alignRight 14 nReject
    ra = alignRight 14 acceptRate
    ro = alignRight 14 optimalRate
    tp = alignRight 20 tuneParam
    mt = alignRight 30 manualAdjustment

-- | For each key @k@, store the number of accepted and rejected proposals.
newtype Acceptance k = Acceptance {fromAcceptance :: M.Map k (Int, Int)}
  deriving (Eq, Read, Show)

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
pushA :: Ord k => k -> Bool -> Acceptance k -> Acceptance k
pushA k True = Acceptance . M.adjust (force . first succ) k . fromAcceptance
pushA k False = Acceptance . M.adjust (force . second succ) k . fromAcceptance
{-# INLINEABLE pushA #-}

-- | Reset acceptance storage.
resetA :: Ord k => Acceptance k -> Acceptance k
resetA = emptyA . M.keys . fromAcceptance

transformKeys :: (Ord k1, Ord k2) => [k1] -> [k2] -> M.Map k1 v -> M.Map k2 v
transformKeys ks1 ks2 m = foldl' insrt M.empty $ zip ks1 ks2
  where
    insrt m' (k1, k2) = M.insert k2 (m M.! k1) m'

-- | Transform keys using the given lists. Keys not provided will not be present
-- in the new 'Acceptance' variable.
transformKeysA :: (Ord k1, Ord k2) => [k1] -> [k2] -> Acceptance k1 -> Acceptance k2
transformKeysA ks1 ks2 = Acceptance . transformKeys ks1 ks2 . fromAcceptance

-- | Acceptance counts and rate for a specific proposal.
--
-- Return 'Nothing' if no proposals have been accepted or rejected (division by
-- zero).
acceptanceRate :: Ord k => k -> Acceptance k -> Maybe (Int, Int, Double)
acceptanceRate k a = case fromAcceptance a M.!? k of
  Just (0, 0) -> Nothing
  Just (as, rs) -> Just (as, rs, fromIntegral as / fromIntegral (as + rs))
  Nothing -> error "acceptanceRate: Key not found in map."

-- | Acceptance rates for all proposals.
--
-- Set rate to 'Nothing' if no proposals have been accepted or rejected
-- (division by zero).
acceptanceRates :: Acceptance k -> M.Map k (Maybe Double)
acceptanceRates =
  M.map
    ( \(as, rs) ->
        if as + rs == 0
          then Nothing
          else Just $ fromIntegral as / fromIntegral (as + rs)
    )
    . fromAcceptance

-- | Header of proposal summaries.
proposalHeader :: BL.ByteString
proposalHeader =
  renderRow
    "Name"
    "Description"
    "Weight"
    "Accepted"
    "Rejected"
    "Rate"
    "Optimal rate"
    "Tuning parameter"
    "Consider manual adjustment"

-- | Horizontal line of proposal summaries.
proposalHLine :: BL.ByteString
proposalHLine = BL.replicate (BL.length proposalHeader) '-'

-- | Proposal summary.
summarizeProposal ::
  PName ->
  PDescription ->
  PWeight ->
  Maybe TuningParameter ->
  PDimension ->
  Maybe (Int, Int, Double) ->
  BL.ByteString
summarizeProposal name description weight tuningParameter dimension ar =
  renderRow
    (BL.pack $ fromPName name)
    (BL.pack $ fromPDescription description)
    weightStr
    nAccept
    nReject
    acceptRate
    optimalRate
    tuneParamStr
    manualAdjustmentStr
  where
    weightStr = BB.toLazyByteString $ BB.intDec $ fromPWeight weight
    nAccept = BB.toLazyByteString $ maybe "" (BB.intDec . (^. _1)) ar
    nReject = BB.toLazyByteString $ maybe "" (BB.intDec . (^. _2)) ar
    acceptRate = BL.fromStrict $ maybe "" (BC.toFixed 2 . (^. _3)) ar
    optimalRate = BL.fromStrict $ BC.toFixed 2 $ getOptimalRate dimension
    tuneParamStr = BL.fromStrict $ maybe "" (BC.toFixed 3) tuningParameter
    checkRate rate
      | rate < rateMin = Just "rate too low"
      | rate > rateMax = Just "rate too high"
      | otherwise = Nothing
    checkTuningParam tp
      | tp <= (1.1 * tuningParameterMin) = Just "tuning parameter too low"
      | tp >= (0.9 * tuningParameterMax) = Just "tuning parameter too high"
      | otherwise = Nothing
    tps = checkTuningParam =<< tuningParameter
    ars = (checkRate . (^. _3)) =<< ar
    manualAdjustmentStr =
      let
       in case (ars, tps) of
            (Nothing, Nothing) -> ""
            (Just s, _) -> s
            (_, Just s) -> s

-- | Summarize the 'Proposal's in the 'Cycle'. Also report acceptance rates.
summarizeCycle :: Acceptance (Proposal a) -> Cycle a -> BL.ByteString
summarizeCycle a c =
  BL.intercalate "\n" $
    [ "Summary of proposal(s) in cycle.",
      nProposalsFullStr,
      describeOrder (ccOrder c),
      proposalHeader,
      proposalHLine
    ]
      ++ [ summarizeProposal
             (pName p)
             (pDescription p)
             (pWeight p)
             (tParam <$> pTuner p)
             (pDimension p)
             (ar p)
           | p <- ps
         ]
      ++ [proposalHLine]
  where
    ps = ccProposals c
    nProposals = getNProposalsPerCycle c
    nProposalsStr = BB.toLazyByteString $ BB.intDec nProposals
    nProposalsFullStr = case nProposals of
      1 -> nProposalsStr <> " proposal is performed per iteration."
      _ -> nProposalsStr <> " proposals are performed per iterations."
    ar m = acceptanceRate m a
