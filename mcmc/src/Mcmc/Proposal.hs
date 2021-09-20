{-# LANGUAGE DerivingVia #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Proposal
-- Description :  Proposals are instruction to move around the state space
-- Copyright   :  (c) Dominik Schrempf 2021
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
    Proposal (..),
    KernelRatio,
    Jacobian,
    JacobianFunction,
    (@~),
    liftProposal,
    liftProposalWith,
    ProposalSimple,
    createProposal,

    -- * Tuners
    Tuner (..),
    Tune (..),
    TuningParameter,
    TuningFunction,
    AuxiliaryTuningParameters,
    AuxiliaryTuningFunction,
    defaultTuningFunctionWith,
    noTuningFunction,
    noAuxiliaryTuningFunction,
    tuningParameterMin,
    tuningParameterMax,
    tuneWithTuningParameters,
    getOptimalRate,

    -- * Output
    proposalHeader,
    summarizeProposal,
  )
where

import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Double.Conversion.ByteString as BC
import Data.Function
import qualified Data.Vector as VB
import qualified Data.Vector.Unboxed as VU
import Lens.Micro
import Mcmc.Acceptance
import Mcmc.Internal.ByteString
import Numeric.Log hiding (sum)
import System.Random.MWC

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
    -- | Dimension of the proposal. The dimension is used to calculate the
    -- optimal acceptance rate, and does not have to be exact.
    prDimension :: PDimension,
    -- | The weight determines how often a 'Proposal' is executed per iteration of
    -- the Markov chain.
    prWeight :: PWeight,
    -- | Simple proposal without name, weight, and tuning information.
    prSimple :: ProposalSimple a,
    -- | Tuning is disabled if set to 'Nothing'.
    prTuner :: Maybe (Tuner a)
  }

instance Eq (Proposal a) where
  m == n = prName m == prName n && prDescription m == prDescription n

instance Ord (Proposal a) where
  compare = compare `on` (\p -> (prDescription p, prName p, prWeight p))

-- | Ratio of the proposal kernels.
--
-- Part of the MHG acceptance ratio.
--
-- See also 'Jacobian'.
--
-- NOTE: Actually the 'Jacobian' should be part of the 'KernelRatio'. However,
-- it is more declarative to have them separate. Like so, we are constantly
-- reminded: Is the Jacobian modifier different from 1.0?
type KernelRatio = Log Double

-- | Absolute value of the determinant of the Jacobian matrix.
--
-- Part of the MHG acceptance ratio.
--
-- See also 'Jacobian'.
type Jacobian = Log Double

-- | Function calculating the 'Jacobian' of a proposal.
type JacobianFunction a = a -> Jacobian

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
-- A function to calculate the Jacobian has to be provided (but see '(@~)').
--
-- For further reference, please see the [example
-- @Pair@](https://github.com/dschrempf/mcmc/blob/master/mcmc-examples/Pair/Pair.hs).
liftProposalWith :: JacobianFunction b -> Lens' b a -> Proposal a -> Proposal b
liftProposalWith jf l (Proposal n r d w s t) =
  Proposal n r d w (liftProposalSimpleWith jf l s) (liftTunerWith jf l <$> t)

-- | Simple proposal without tuning information.
--
-- Instruction about randomly moving from the current state to a new state,
-- given some source of randomness.
--
-- In order to calculate the Metropolis-Hastings-Green ratio, we need to know
-- the ratio of the backward to forward kernels (the 'KernelRatio' or the
-- probability masses or probability densities) and the 'Jacobian'.
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
type ProposalSimple a = a -> GenIO -> IO (a, KernelRatio, Jacobian)

-- Lift a simple proposal from one data type to another.
liftProposalSimpleWith :: JacobianFunction b -> Lens' b a -> ProposalSimple a -> ProposalSimple b
liftProposalSimpleWith jf l s = s'
  where
    s' y g = do
      (x', r, j) <- s (y ^. l) g
      let y' = set l x' y
          jxy = jf y
          jyx = jf y'
          j' = j * jyx / jxy
      return (y', r, j')

-- | Create a proposal with a single tuning parameter.
--
-- Proposals with arbitrary tuning parameters have to be created manually. See
-- 'Tuner' for more information, and 'Mcmc.Proposal.Hamiltonian' for an example.
createProposal ::
  -- | Description of the proposal type and parameters.
  PDescription ->
  -- | Function creating a simple proposal for a given tuning parameter.
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
createProposal r f d n w Tune =
  Proposal n r d w (f 1.0) (Just tuner)
  where
    fT = defaultTuningFunctionWith d
    fTs = noAuxiliaryTuningFunction
    g t _ = Right $ f t
    tuner = Tuner 1.0 fT VU.empty fTs g
createProposal r f d n w NoTune =
  Proposal n r d w (f 1.0) Nothing

-- | Required information to tune 'Proposal's.
data Tuner a = Tuner
  { tTuningParameter :: TuningParameter,
    tTuningFunction :: TuningFunction,
    tAuxiliaryTuningParameters :: AuxiliaryTuningParameters,
    tAuxiliaryTuningFunction :: AuxiliaryTuningFunction a,
    -- | Given the tuning parameter, and the auxiliary tuning parameters, get
    -- the tuned simple proposal.
    --
    -- Should return 'Left' if the vector of auxiliary tuning parameters is
    -- invalid.
    tGetSimpleProposal ::
      TuningParameter ->
      AuxiliaryTuningParameters ->
      Either String (ProposalSimple a)
  }

-- Lift tuner from one data type to another.
liftTunerWith :: JacobianFunction b -> Lens' b a -> Tuner a -> Tuner b
liftTunerWith jf l (Tuner p fP ps fPs g) = Tuner p fP ps fPs' g'
  where
    fPs' = fPs . VB.map (^. l)
    g' x xs = liftProposalSimpleWith jf l <$> g x xs

-- | Tune proposal?
data Tune = Tune | NoTune
  deriving (Show, Eq)

-- | Tuning parameter.
--
--  The larger the tuning parameter, the larger the proposal and the lower the
-- expected acceptance rate; and vice versa.
type TuningParameter = Double

-- | Compute new tuning parameter from a given acceptance rate and the old
-- tuning parameter.
type TuningFunction = AcceptanceRate -> TuningParameter -> TuningParameter

-- | Auxiliary tuning parameters; vector may be empty.
--
-- Auxiliary tuning parameters are not shown in proposal summaries.
type AuxiliaryTuningParameters = VU.Vector TuningParameter

-- | Compute new auxiliary tuning parameters from a given trace and the old
-- auxiliary tuning parameters.
type AuxiliaryTuningFunction a = VB.Vector a -> AuxiliaryTuningParameters -> AuxiliaryTuningParameters

-- | Default tuning function.
--
-- Subject to change.
defaultTuningFunctionWith ::
  -- Optimal acceptance rate.
  PDimension ->
  TuningFunction
defaultTuningFunctionWith d r t = let rO = getOptimalRate d in exp (2 * (r - rO)) * t

-- | Do not tune.
--
-- Useful if auxiliary tuning parameters are tuned, but not the main tuning
-- parameter.
noTuningFunction :: TuningFunction
noTuningFunction _ = id

-- | Do not tune auxiliary parameters.
noAuxiliaryTuningFunction :: AuxiliaryTuningFunction a
noAuxiliaryTuningFunction _ = id

-- IDEA: Per proposal type tuning parameter boundaries. For example, a sliding
-- proposal with a large tuning parameter is not a problem. But then, if the
-- tuning parameters are very different from one, a different base proposal
-- should be chosen.

-- | Minimal tuning parameter; subject to change.
tuningParameterMin :: TuningParameter
tuningParameterMin = 1e-5

-- | Maximal tuning parameter; subject to change.
tuningParameterMax :: TuningParameter
tuningParameterMax = 1e3

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
  Just (Tuner _ fT _ fTs g) ->
    -- Ensure that the tuning parameter is strictly positive and well bounded.
    let t' = max tuningParameterMin t
        t'' = min tuningParameterMax t'
        psE = g t'' ts
     in case psE of
          Left err -> Left $ "tune: " <> err
          Right ps -> Right $ p {prSimple = ps, prTuner = Just $ Tuner t'' fT ts fTs g}

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
    tuneParamStr = BL.fromStrict $ maybe "" (BC.toFixed 4) tuningParameter
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
