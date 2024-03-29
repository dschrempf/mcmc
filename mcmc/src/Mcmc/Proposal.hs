{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DerivingVia #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal
-- Description :  Proposals are instruction to move around the state space
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 13:42:53 2020.
module Mcmc.Proposal
  ( -- * Proposals
    PName (..),
    PDescription (..),
    PWeight (fromPWeight),
    pWeight,
    PDimension (..),
    PSpeed (..),
    Proposal (..),
    KernelRatio,
    PResult (..),
    Jacobian,
    JacobianFunction,
    PFunction,
    createProposal,

    -- * Tuners
    Tuner (..),
    Tune (..),
    TuningParameter,
    TuningType (..),
    TuningFunction,
    AuxiliaryTuningParameters,
    tuningFunction,
    tuningParameterMin,
    tuningParameterMax,
    tuneWithTuningParameters,
    getOptimalRate,

    -- * Lift proposals
    (@~),
    liftProposal,
    liftProposalWith,

    -- * Output
    proposalHeader,
    summarizeProposal,
  )
where

import Data.Bifunctor
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Function
import qualified Data.Vector as VB
import qualified Data.Vector.Unboxed as VU
import Lens.Micro
import Mcmc.Acceptance
import Mcmc.Internal.ByteString
import Mcmc.Jacobian
import Numeric.Log hiding (sum)
import System.Random.Stateful

-- | Proposal name.
newtype PName = PName {fromPName :: String}
  deriving (Show, Eq, Ord)
  deriving (Monoid, Semigroup) via String

-- | Proposal description.
newtype PDescription = PDescription {fromPDescription :: String}
  deriving (Show, Eq, Ord)

-- | The positive weight determines how often a 'Proposal' is executed per
-- iteration of the Markov chain. Abstract data type; for construction, see
-- 'pWeight'.
newtype PWeight = PWeight {fromPWeight :: Int}
  deriving (Show, Eq, Ord)

-- | Check if the weight is positive.
--
-- Call 'error' if weight is zero or negative.
pWeight :: Int -> PWeight
pWeight n
  | n <= 0 = error $ "pWeight: Proposal weight is zero or negative: " <> show n <> "."
  | otherwise = PWeight n

-- | Proposal dimension.
--
-- The number of affected, independent parameters.
--
-- The dimension is used to calculate the optimal acceptance rate, and does not
-- have to be exact.
--
-- Usually, the optimal acceptance rate of low dimensional proposals is higher
-- than for high dimensional ones. However, this is not always true (see below).
--
-- Further, optimal acceptance rates are still subject to controversies. To my
-- knowledge, research has focused on random walk proposals with multivariate
-- normal distributions of dimension @d@. In this case, the following acceptance
-- rates are desired:
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
-- rescaling all other variables. What is the dimension of this proposal? Here,
-- the dimension is set to 2. The reason is that if the dimension of the simplex
-- is 2, two variables are changed. If the dimension of the simplex is high, one
-- variable is changed substantially, while all others are changed marginally.
--
-- Further, if a proposal changes a number of variables in the same way (and not
-- independently like in a random walk proposal), the dimension of the proposal
-- is set to the number of variables changed.
--
-- Moreover, proposals of unknown dimension are assumed to have high dimension,
-- and the optimal acceptance rate 0.234 is used.
--
-- Finally, special proposals may have completely different desired acceptance
-- rates. For example. the Hamiltonian Monte Carlo proposal (see
-- Mcmc.Proposal.Hamiltonian.hmc) has a desired acceptance rate of 0.65.
-- Specific acceptance rates can be set with 'PSpecial'.
data PDimension
  = PDimension Int
  | PDimensionUnknown
  | -- | Provide dimension ('Int') and desired acceptance rate ('Double').
    PSpecial Int Double

-- | Rough indication whether a proposal is fast or slow.
--
-- Useful during burn in. Slow proposals are not executed during fast auto
-- tuning periods.
--
-- See 'Mcmc.Settings.BurnInSettings'.
data PSpeed = PFast | PSlow
  deriving (Eq)

-- | A 'Proposal' is an instruction about how the Markov chain will traverse the
-- state space @a@. Essentially, it is a probability mass or probability density
-- conditioned on the current state (i.e., a Markov kernel).
--
-- A 'Proposal' may be tuneable in that it contains information about how to
-- enlarge or shrink the proposal size to decrease or increase the acceptance
-- rate.
--
-- Predefined proposals are provided. To create custom proposals, one may use
-- the convenience function 'createProposal'.
data Proposal a = Proposal
  { -- | Name of the affected variable.
    prName :: PName,
    -- | Description of the proposal type and parameters.
    prDescription :: PDescription,
    prSpeed :: PSpeed,
    prDimension :: PDimension,
    -- | The weight determines how often a 'Proposal' is executed per iteration of
    -- the Markov chain.
    prWeight :: PWeight,
    -- | Simple proposal function without name, weight, and tuning information.
    prFunction :: PFunction a,
    -- | Tuning is disabled if set to 'Nothing'.
    prTuner :: Maybe (Tuner a)
  }

instance Eq (Proposal a) where
  m == n = prName m == prName n && prDescription m == prDescription n

instance Ord (Proposal a) where
  compare = compare `on` (\p -> (prDescription p, prName p, prWeight p))

-- | Ratio of the proposal kernels.
--
-- For unbiased, volume preserving proposals, the values is 1.0.
--
-- For biased proposals, the kernel ratio is qYX / qXY, where qAB is the
-- probability density to move from A to B.
type KernelRatio = Log Double

-- | Proposal result.
data PResult a
  = -- | Accept the new value regardless of the prior, likelihood or Jacobian.
    ForceAccept !a
  | -- | Reject the proposal regardless of the prior, likelihood or Jacobian.
    ForceReject
  | -- | Propose a new value.
    --
    -- In order to calculate the Metropolis-Hastings-Green ratio, we need to know
    -- the ratio of the backward to forward kernels (the 'KernelRatio' or the
    -- probability masses or probability densities) and the 'Jacobian'.
    --
    -- Note: The 'Jacobian' should be part of the 'KernelRatio'. However, it is
    -- more declarative to have them separate. Like so, we are constantly
    -- reminded: Is the Jacobian modifier different from 1.0?
    Propose !a !KernelRatio !Jacobian
  deriving (Show, Eq)

-- | Simple proposal function without tuning information.
--
-- Instruction about randomly moving from the current state to a new state,
-- given some source of randomness.
--
-- Maybe report acceptance rates internal to the proposal (e.g., used by
-- proposals based on Hamiltonian dynamics).
type PFunction a = a -> IOGenM StdGen -> IO (PResult a, Maybe AcceptanceRates)

-- | Create a proposal with a single tuning parameter.
--
-- Proposals with auxiliary tuning parameters have to be created manually. See
-- 'Tuner' for more information, and 'Mcmc.Proposal.Hamiltonian' for an example.
createProposal ::
  -- | Description of the proposal type and parameters.
  PDescription ->
  -- | Function creating a simple proposal function for a given tuning parameter.
  (TuningParameter -> PFunction a) ->
  -- | Speed.
  PSpeed ->
  -- | Dimension.
  PDimension ->
  -- | Name.
  PName ->
  -- | Weight.
  PWeight ->
  -- | Activate tuning?
  Tune ->
  Proposal a
createProposal r f s d n w Tune =
  Proposal n r s d w (f 1.0) (Just tuner)
  where
    fT = tuningFunction
    g t _ = Right $ f t
    tuner = Tuner 1.0 VU.empty False False fT g
createProposal r f s d n w NoTune =
  Proposal n r s d w (f 1.0) Nothing

-- | Required information to tune 'Proposal's.
data Tuner a = Tuner
  { tTuningParameter :: TuningParameter,
    tAuxiliaryTuningParameters :: AuxiliaryTuningParameters,
    -- | Does the tuner require the trace over the last tuning period?
    tRequireTrace :: Bool,
    -- | Can the tuner be used for intermediate tuning (see 'TuningType')?
    tSuitableForIntermediateTuning :: Bool,
    tTuningFunction :: TuningFunction a,
    -- | Given the tuning parameter, and the auxiliary tuning parameters, get
    -- the tuned propose function.
    --
    -- Should return 'Left' if the vector of auxiliary tuning parameters is
    -- invalid.
    tGetPFunction ::
      TuningParameter ->
      AuxiliaryTuningParameters ->
      Either String (PFunction a)
  }

-- | Tune proposal?
data Tune = Tune | NoTune
  deriving (Show, Eq)

-- | Tuning parameter.
--
--  The larger the tuning parameter, the larger the proposal and the lower the
-- expected acceptance rate; and vice versa.
type TuningParameter = Double

-- | Tuning type. To distinguish between fast and slow proposals, see
-- 'Mcmc.Cycle.IterationMode'.
data TuningType
  = -- | Normal tuning step with fast proposals only.
    NormalTuningFastProposalsOnly
  | -- | Intermediate tuning step executed after each iteration with fast
    -- proposals only. Only suitable for proposals which can calculate expected
    -- acceptance rates.
    IntermediateTuningFastProposalsOnly
  | -- | The last tuning step with fast proposals only may be special.
    LastTuningFastProposalsOnly
  | -- | Normal tuning step of all proposals.
    NormalTuningAllProposals
  | -- | Intermediate tuning step of all proposals.
    IntermediateTuningAllProposals
  | -- | The last tuning step with all proposals.
    LastTuningAllProposals
  deriving (Eq)

-- | Compute new tuning parameters.
type TuningFunction a =
  TuningType ->
  PDimension ->
  -- | Acceptance rate of last tuning period. May not always be available
  -- because proposals may be skipped.
  Maybe AcceptanceRate ->
  -- | Trace of last tuning period. Not available for intermediate tuning' steps
  -- (see 'TuningType'), and only available for other tuning types when
  -- requested by proposal.
  Maybe (VB.Vector a) ->
  (TuningParameter, AuxiliaryTuningParameters) ->
  (TuningParameter, AuxiliaryTuningParameters)

-- | Auxiliary tuning parameters.
--
-- Auxiliary tuning parameters are not shown in proposal summaries.
--
-- Vector may be empty.
type AuxiliaryTuningParameters = VU.Vector TuningParameter

tuningFunctionSimple :: PDimension -> AcceptanceRate -> TuningParameter -> TuningParameter
tuningFunctionSimple d r t = let rO = getOptimalRate d in exp (2 * (r - rO)) * t

-- | Default tuning function.
--
-- The default tuning function only uses the actual acceptance rate. In
-- particular, it does not handle auxiliary tuning parameters, ignores
-- intermediate tuning steps, and ignores the actual samples attained during the
-- last tuning period.
tuningFunction :: TuningFunction a
tuningFunction IntermediateTuningFastProposalsOnly _ _ _ t = t
tuningFunction IntermediateTuningAllProposals _ _ _ t = t
tuningFunction _ _ Nothing _ t = t
tuningFunction _ d (Just r) _ (!t, !ts) = first (tuningFunctionSimple d r) (t, ts)

-- IDEA: Per proposal type tuning parameter boundaries. For example, a sliding
-- proposal with a large tuning parameter is not a problem. But then, if the
-- tuning parameters are very different from one, a different base proposal
-- should be chosen.

-- | Minimal tuning parameter; subject to change.
tuningParameterMin :: TuningParameter
tuningParameterMin = 1e-6

-- | Maximal tuning parameter; subject to change.
tuningParameterMax :: TuningParameter
tuningParameterMax = 5e3

-- | Tune a 'Proposal'.
--
-- The size of the proposal is proportional to the tuning parameter which has
-- positive lower and upper boundaries of 'tuningParameterMin' and
-- 'tuningParameterMax', respectively.
--
-- Auxiliary tuning parameters may also be used by the 'Tuner' of the proposal.
--
-- Return 'Left' if:
--
-- - the 'Proposal' is not tuneable;
--
-- - the auxiliary tuning parameters are invalid.
--
-- Used by 'Mcmc.Chain.Save.fromSavedChain'.
tuneWithTuningParameters ::
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Proposal a ->
  Either String (Proposal a)
tuneWithTuningParameters t ts p = case prTuner p of
  Nothing -> Left "tuneWithTuningParameters: Proposal is not tunable."
  Just (Tuner _ _ reqTr inTn fT g) ->
    -- Ensure that the tuning parameter is strictly positive and well bounded.
    let t' = max tuningParameterMin t
        t'' = min tuningParameterMax t'
        psE = g t'' ts
     in case psE of
          Left err -> Left $ "tune: " <> err
          Right ps -> Right $ p {prFunction = ps, prTuner = Just $ Tuner t'' ts reqTr inTn fT g}

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
getOptimalRate (PSpecial _ r) = r

-- | Lift a proposal from one data type to another.
--
-- Assume the Jacobian is 1.0.
--
-- For example:
--
-- @
-- scaleFirstEntryOfTuple = _1 @~ scale
-- @
--
-- See also 'liftProposalWith'.
infixl 7 @~

(@~) :: Lens' b a -> Proposal a -> Proposal b
(@~) = liftProposal

-- | See '(@~)'.
liftProposal :: Lens' b a -> Proposal a -> Proposal b
liftProposal = liftProposalWith (const 1.0)

-- | Lift a proposal from one data type to another.
--
-- A function to calculate the Jacobian has to be provided. If the Jabobian is
-- 1.0, use 'liftProposal'.
--
-- For further reference, please see the [example
-- @Pair@](https://github.com/dschrempf/mcmc/blob/master/mcmc-examples/Pair/Pair.hs).
liftProposalWith :: JacobianFunction b -> Lens' b a -> Proposal a -> Proposal b
liftProposalWith jf l (Proposal n r d p w s t) =
  Proposal n r d p w (liftPFunctionWith jf l s) (liftTunerWith jf l <$> t)

-- Lift a proposal function from one data type to another.
liftPFunctionWith :: JacobianFunction b -> Lens' b a -> PFunction a -> PFunction b
liftPFunctionWith jf l s = s'
  where
    s' y g = do
      (pr, ac) <- s (y ^. l) g
      let pr' = case pr of
            ForceAccept x' -> ForceAccept $ set l x' y
            ForceReject -> ForceReject
            Propose x' r j ->
              let y' = set l x' y
                  jxy = jf y
                  jyx = jf y'
                  j' = j * jyx / jxy
               in Propose y' r j'
      pure (pr', ac)

-- Lift tuner from one data type to another.
liftTunerWith :: JacobianFunction b -> Lens' b a -> Tuner a -> Tuner b
liftTunerWith jf l (Tuner p ps reqTr inTn fP g) = Tuner p ps reqTr inTn fP' g'
  where
    fP' b d r = fP b d r . fmap (VB.map (^. l))
    g' x xs = liftPFunctionWith jf l <$> g x xs

-- Warn if acceptance rate is lower.
rateMin :: Double
rateMin = 0.1

-- Warn if acceptance rate is larger.
rateMax :: Double
rateMax = 0.9

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
renderRow name ptype weight nAccept nReject acceptRateActual optimalRate tuneParam manualAdjustment =
  nm <> pt <> wt <> na <> nr <> ra <> ro <> tp <> mt
  where
    nm = alignLeft 30 name
    pt = alignLeft 50 ptype
    wt = alignRight 8 weight
    na = alignRight 14 nAccept
    nr = alignRight 14 nReject
    ra = alignRight 14 acceptRateActual
    ro = alignRight 14 optimalRate
    tp = alignRight 20 tuneParam
    mt = alignRight 30 manualAdjustment

-- | Header of proposal summaries.
proposalHeader :: BL.ByteString
proposalHeader =
  renderRow
    "Name"
    "Description"
    "Weight"
    "Accepted"
    "Rejected"
    "Actual rate"
    "Optimal rate"
    "Tuning parameter"
    "Consider manual adjustment"

-- | Proposal summary.
summarizeProposal ::
  PName ->
  PDescription ->
  PWeight ->
  Maybe TuningParameter ->
  PDimension ->
  (Int, Int, Maybe Double, Maybe Double) ->
  BL.ByteString
summarizeProposal name description weight tuningParameter dimension ar =
  renderRow
    (BL.pack $ fromPName name)
    (BL.pack $ fromPDescription description)
    weightStr
    nAccept
    nReject
    acceptRateActual
    -- acceptRateExpected
    optimalRate
    tuneParamStr
    manualAdjustmentStr
  where
    fN n = BB.formatDouble (BB.standard n)
    weightStr = BB.toLazyByteString $ BB.intDec $ fromPWeight weight
    nAccept = BB.toLazyByteString $ BB.intDec $ ar ^. _1
    nReject = BB.toLazyByteString $ BB.intDec $ ar ^. _2
    acceptRateActual = BB.toLazyByteString $ maybe "" (fN 2) (ar ^. _3)
    optimalRate = BB.toLazyByteString $ fN 2 $ getOptimalRate dimension
    tuneParamStr = BB.toLazyByteString $ maybe "" (fN 6) tuningParameter
    checkRate rate
      | rate < rateMin = Just "rate too low"
      | rate > rateMax = Just "rate too high"
      | otherwise = Nothing
    checkTuningParam tp
      | tp <= (1.1 * tuningParameterMin) = Just "tuning parameter too low"
      | tp >= (0.9 * tuningParameterMax) = Just "tuning parameter too high"
      | otherwise = Nothing
    tps = checkTuningParam =<< tuningParameter
    -- Use actual acceptance rate.
    ars = checkRate =<< (ar ^. _3)
    manualAdjustmentStr =
      let
       in case (ars, tps) of
            (Nothing, Nothing) -> ""
            (Just s, _) -> s
            (_, Just s) -> s
