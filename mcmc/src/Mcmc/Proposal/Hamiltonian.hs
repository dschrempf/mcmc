{-# LANGUAGE TupleSections #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      :  Mcmc.Proposal.Hamiltonian
-- Description :  Hamiltonian Monte Carlo proposal
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Mon Jul  5 12:59:42 2021.
--
-- The Hamiltonian Monte Carlo (HMC) proposal.
--
-- For references, see:
--
-- - [1] Chapter 5 of Handbook of Monte Carlo: Neal, R. M., MCMC Using
--   Hamiltonian Dynamics, In S. Brooks, A. Gelman, G. Jones, & X. Meng (Eds.),
--   Handbook of Markov Chain Monte Carlo (2011), CRC press.
--
-- - [2] Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B., Bayesian data
--   analysis (2014), CRC Press.
--
-- - [3] Review by Betancourt and notes: Betancourt, M., A conceptual
--   introduction to Hamiltonian Monte Carlo, arXiv, 1701–02434 (2017).
--
-- NOTE on implementation:
--
-- - The implementation assumes the existence of the gradient. Like so, the user
--   can use automatic or manual differentiation, depending on the problem at
--   hand.
--
-- - The state needs to be list like or 'Traversable' so that the structure of
--   the state space is available. A 'Traversable' constraint on the data type
--   is nice because it is more general than, for example, a list, and specific
--   data structures can be used.
--
-- - The state needs to have a zip-like 'Applicative' instance so that
-- - matrix/vector operations can be performed.

module Mcmc.Proposal.Hamiltonian
  ( Gradient,
    Masses,
    LeapfrogTrajectoryLength,
    LeapfrogScalingFactor,
    HTune (..),
    HSettings (..),
    hamiltonian,
  )
where

import Data.Foldable
import qualified Data.Matrix as M
import Data.Maybe
import Data.Traversable
import qualified Data.Vector as VB
import Mcmc.Prior
import Mcmc.Proposal
import Numeric.Log
import Statistics.Distribution
import Statistics.Distribution.Normal
import qualified Statistics.Function as S
import qualified Statistics.Sample as S
import System.Random.MWC

-- TODO: At the moment, the HMC proposal is agnostic of the prior and
-- likelihood, that is, the posterior function. This means, that it cannot know
-- when it reaches a point with zero posterior probability. This also affects
-- restricted or constrained parameters. See Gelman p. 303.

-- TODO: No-U-turn sampler.

-- TODO: Riemannian adaptation.

-- | Gradient of the log posterior function.
type Gradient f = f Double -> f Double

-- | Function validating the state.
--
-- Useful if parameters are constrained.
type Validate f = f Double -> Bool

-- | Masses of parameters.
--
-- NOTE: Full specification of a mass matrix including off-diagonal elements is
-- not supported.
--
-- NOTE: Parameters without masses ('Nothing') are not changed by the
-- Hamiltonian proposal.
--
-- The masses roughly describe how reluctant the particle moves through the
-- state space. If a parameter has higher mass, the momentum in this direction
-- will be changed less by the provided gradient, than when the same parameter
-- has lower mass.
--
-- The proposal is more efficient if masses are assigned according to the
-- inverse (co)-variance structure of the posterior function. That is,
-- parameters changing on larger scales should have lower masses than parameters
-- changing on lower scales. In particular, and for a diagonal mass matrix, the
-- optimal masses are the inverted variances of the parameters distributed
-- according to the posterior function.
--
-- Of course, the scales of the parameters of the posterior function are usually
-- unknown. Often, it is sufficient to
--
-- - set the masses to identical values roughly scaled with the inverted
--   estimated average variance of the posterior function; or even to
--
-- - set all masses to 1.0, and trust the tuning algorithm (see
--   'HTuneMassesAndLeapfrog') to find the correct values.
type Masses f = f (Maybe Double)

-- | Mean leapfrog trajectory length \(L\).
--
-- Number of leapfrog steps per proposal.
--
-- To avoid problems with ergodicity, the actual number of leapfrog steps is
-- sampled proposal from a discrete uniform distribution over the interval
-- \([\text{floor}(0.8L),\text{ceiling}(1.2L)]\).
--
-- For a discussion of ergodicity and reasons why randomization is important,
-- see [1] p. 15; also mentioned in [2] p. 304.
--
-- NOTE: To avoid errors, the left bound has an additional hard minimum of 1,
-- and the right bound is required to be larger equal than the left bound.
--
-- Usually set to 10, but larger values may be desirable.
type LeapfrogTrajectoryLength = Int

-- | Mean of leapfrog scaling factor \(\epsilon\).
--
-- Determines the size of each leapfrog step.
--
-- To avoid problems with ergodicity, the actual leapfrog scaling factor is
-- sampled per proposal from a continuous uniform distribution over the interval
-- \((0.8\epsilon,1.2\epsilon]\).
--
-- For a discussion of ergodicity and reasons why randomization is important,
-- see [1] p. 15; also mentioned in [2] p. 304.
--
-- Usually set such that \( L \epsilon = 1.0 \), but smaller values may be
-- required if acceptance rates are low.
type LeapfrogScalingFactor = Double

-- Target state containing parameters.
type Positions f = f Double

-- Momenta of the parameters.
type Momenta f = f (Maybe Double)

-- | Tuning settings.
--
--
-- Tuning of leapfrog parameters:
--
-- We expect that the larger the leapfrog step size the larger the proposal step
-- size and the lower the acceptance ratio. Consequently, if the acceptance rate
-- is too low, the leapfrog step size is decreased and vice versa. Further, the
-- leapfrog trajectory length is scaled such that the product of the leapfrog
-- step size and trajectory length stays constant.
--
-- Tuning of masses:
--
-- The variances of all parameters of the posterior distribution obtained over
-- the last auto tuning interval is calculated and the masses are amended using
-- the old masses and the inverted variances. If, for a specific coordinate, the
-- sample size is too low, or if the calculated variance is out of predefined
-- bounds, the mass of the affected position is not changed.
data HTune
  = -- | Tune masses and leapfrog parameters.
    HTuneMassesAndLeapfrog
  | -- | Tune leapfrog parameters only.
    HTuneLeapfrogOnly
  | -- | Do not tune at all.
    HNoTune
  deriving (Eq, Show)

-- | Specifications for Hamilton Monte Carlo proposal.
data HSettings f = HSettings
  { hGradient :: Gradient f,
    hMaybeValidate :: Maybe (Validate f),
    hMasses :: Masses f,
    hLeapfrogTrajectoryLength :: LeapfrogTrajectoryLength,
    hLeapfrogScalingFactor :: LeapfrogScalingFactor,
    hTune :: HTune
  }

checkHSettings :: Foldable f => HSettings f -> Maybe String
checkHSettings (HSettings _ _ masses l eps _)
  | any f masses = Just "checkHSettings: One or more masses are zero or negative."
  | l < 1 = Just "checkHSettings: Leapfrog trajectory length is zero or negative."
  | eps <= 0 = Just "checkHSettings: Leapfrog scaling factor is zero or negative."
  | otherwise = Nothing
  where
    f (Just m) = m <= 0
    f Nothing = False

generateMomenta ::
  Traversable f =>
  Masses f ->
  GenIO ->
  IO (Momenta f)
generateMomenta masses gen = traverse (generateWith gen) masses
  where
    generateWith g (Just m) = let d = normalDistr 0 (sqrt m) in Just <$> genContVar d g
    generateWith _ Nothing = pure Nothing

priorMomenta ::
  (Applicative f, Foldable f) =>
  Masses f ->
  Momenta f ->
  Prior
priorMomenta masses phi = foldl' (*) 1.0 $ f <$> masses <*> phi
  where
    f (Just m) (Just p) = let d = normalDistr 0 (sqrt m) in Exp $ logDensity d p
    f Nothing Nothing = 1.0
    f _ _ = error "priorMomenta: Got just a mass and no momentum or the other way around."

leapfrog ::
  Applicative f =>
  Gradient f ->
  Maybe (Validate f) ->
  Masses f ->
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  Positions f ->
  Momenta f ->
  -- Maybe (Positions', Momenta').
  Maybe (Positions f, Momenta f)
leapfrog grad mVal masses l eps theta phi = do
  let -- The first half step of the momenta.
      phiHalf = leapfrogStepMomenta 0.5 eps grad theta phi
  -- L-1 full steps. This gives the positions theta_{L-1}, and the momenta
  -- phi_{L-1/2}.
  (thetaLM1, phiLM1Half) <- go (l - 1) (Just (theta, phiHalf))
  -- The last full step of the positions.
  thetaL <- valF $ leapfrogStepPositions eps masses thetaLM1 phiLM1Half
  let -- The last half step of the momenta.
      phiL = leapfrogStepMomenta 0.5 eps grad thetaL phiLM1Half
  return (thetaL, phiL)
  where
    valF x = case mVal of
      Nothing -> Just x
      Just f -> if f x then Just x else Nothing
    go _ Nothing = Nothing
    go 0 (Just (t, p)) = Just (t, p)
    go n (Just (t, p)) =
      let t' = leapfrogStepPositions eps masses t p
          p' = leapfrogStepMomenta 1.0 eps grad t' p
          r = (,p') <$> valF t'
       in go (n - 1) r

leapfrogStepMomenta ::
  Applicative f =>
  -- Size of step (half or full step).
  Double ->
  LeapfrogScalingFactor ->
  Gradient f ->
  -- Current positions.
  Positions f ->
  -- Current momenta.
  Momenta f ->
  -- New momenta.
  Momenta f
leapfrogStepMomenta xi eps grad theta phi = phi <+. ((xi * eps) .* grad theta)
  where
    (<+.) :: Applicative f => f (Maybe Double) -> f Double -> f (Maybe Double)
    (<+.) xs ys = f <$> xs <*> ys
    f Nothing _ = Nothing
    f (Just x) y = Just $ x + y

leapfrogStepPositions ::
  Applicative f =>
  LeapfrogScalingFactor ->
  Masses f ->
  -- Current positions.
  Positions f ->
  -- Current momenta.
  Momenta f ->
  Positions f
-- The arguments are flipped to encounter the maybe momentum.
leapfrogStepPositions eps masses theta phi = theta <+. (mScaledReversed .*> phi)
  where
    (<+.) :: Applicative f => f Double -> f (Maybe Double) -> f Double
    (<+.) xs ys = f <$> xs <*> ys
    f x Nothing = x
    f x (Just y) = x + y
    mScaledReversed = (fmap . fmap) ((* eps) . (** (-1))) masses
    (.*>) :: Applicative f => f (Maybe Double) -> f (Maybe Double) -> f (Maybe Double)
    (.*>) xs ys = g <$> xs <*> ys
    g (Just x) (Just y) = Just $ x * y
    g Nothing Nothing = Nothing
    g _ _ = error "leapfrogStepPositions: Got just a mass and no momentum or the other way around."

-- Scalar-vector multiplication.
(.*) :: Applicative f => Double -> f Double -> f Double
(.*) x ys = (* x) <$> ys

-- NOTE: Fixed parameters without mass have a tuning parameter of NaN.
massesToTuningParameters :: Foldable f => Masses f -> AuxiliaryTuningParameters
massesToTuningParameters = VB.fromList . map (fromMaybe nan) . toList
  where
    nan = 0 / 0

-- We need the structure in order to fill it with the given parameters.
tuningParametersToMasses ::
  Traversable f =>
  AuxiliaryTuningParameters ->
  Masses f ->
  Either String (Masses f)
tuningParametersToMasses xs ms =
  if null xs'
    then sequenceA msE
    else Left "tuningParametersToMasses: Too many values."
  where
    (xs', msE) = mapAccumL setValue (VB.toList xs) ms
    setValue [] _ = ([], Left "tuningParametersToMasses: Too few values.")
    -- NOTE: Recover fixed parameters and unset their mass.
    setValue (y : ys) _ = let y' = if isNaN y then Nothing else Just y in (ys, Right y')

hTuningParametersToSettings ::
  Traversable f =>
  TuningParameter ->
  AuxiliaryTuningParameters ->
  HSettings f ->
  Either String (HSettings f)
hTuningParametersToSettings t ts (HSettings g v m l e tn) =
  if tn == HTuneMassesAndLeapfrog
    then case tuningParametersToMasses ts m of
      Left err -> Left err
      Right m' -> Right $ HSettings g v m' lTuned eTuned tn
    else Right $ HSettings g v m lTuned eTuned tn
  where
    -- The larger epsilon, the larger the proposal step size and the lower the
    -- expected acceptance ratio.
    --
    -- Further, we roughly keep \( L * \epsilon = 1.0 \). The equation is not
    -- correct, because we pull L closer to the original value to keep the
    -- runtime somewhat acceptable.
    lTuned = ceiling $ fromIntegral l / (t ** 0.9) :: Int
    eTuned = t * e

hamiltonianSimpleWithTuningParameters ::
  (Applicative f, Traversable f) =>
  HSettings f ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (ProposalSimple (Positions f))
hamiltonianSimpleWithTuningParameters s t ts = case hTuningParametersToSettings t ts s of
  Left err -> Left err
  Right s' -> Right $ hamiltonianSimple s'

hamiltonianSimple ::
  (Applicative f, Traversable f) =>
  HSettings f ->
  ProposalSimple (Positions f)
hamiltonianSimple (HSettings gradient mVal masses l e _) theta g = do
  phi <- generateMomenta masses g
  lRan <- uniformR (lL, lR) g
  eRan <- uniformR (eL, eR) g
  case leapfrog gradient mVal masses lRan eRan theta phi of
    Nothing -> return (theta, 0.0, 1.0)
    Just (theta', phi') ->
      let prPhi = priorMomenta masses phi
          -- NOTE: Neal page 12: In order for the proposal to be in detailed
          -- balance, the momenta have to be negated before proposing the new value.
          -- This is not required here since the prior involves normal distributions
          -- centered around 0. However, if the multivariate normal distribution is
          -- used, it makes a difference.
          prPhi' = priorMomenta masses phi'
          kernelR = prPhi' / prPhi
       in return (theta', kernelR, 1.0)
  where
    lL = maximum [1 :: Int, floor $ (0.8 :: Double) * fromIntegral l]
    lR = maximum [lL, ceiling $ (1.2 :: Double) * fromIntegral l]
    eL = 0.8 * e
    eR = 1.2 * e

minVariance :: Double
minVariance = 1e-6

maxVariance :: Double
maxVariance = 1e6

minSamples :: Int
minSamples = 60

computeAuxiliaryTuningParameters ::
  Foldable f =>
  VB.Vector (Positions f) ->
  AuxiliaryTuningParameters ->
  AuxiliaryTuningParameters
computeAuxiliaryTuningParameters xss ts =
  VB.zipWith (\t -> rescueWith t . calcSamplesAndVariance) ts xssT
  where
    -- TODO: Improve matrix transposition.
    xssT = VB.fromList $ M.toColumns $ M.fromLists $ VB.toList $ VB.map toList xss
    calcSamplesAndVariance xs = (VB.length $ VB.uniq $ S.gsort xs, S.variance xs)
    rescueWith t (sampleSize, var) =
      if var < minVariance || maxVariance < var || sampleSize < minSamples
        then -- then traceShow ("Rescue with " <> show t) t
          t
        else
          let t' = sqrt (t * recip var)
           in -- in traceShow ("Old mass " <> show t <> " new mass " <> show t') t'
              t'

-- | Hamiltonian Monte Carlo proposal.
--
-- The 'Applicative' and 'Traversable' instances are used for element-wise
-- operations.
--
-- Assume a zip-like 'Applicative' instance so that cardinality remains
-- constant.
--
-- NOTE: The desired acceptance rate is 0.65, although the dimension of the
-- proposal is high.
--
-- NOTE: The speed of this proposal can change drastically when tuned because
-- the leapfrog trajectory length is changed.
hamiltonian ::
  (Applicative f, Traversable f) =>
  -- | The sample state is used to calculate the dimension of the proposal.
  f Double ->
  HSettings f ->
  PName ->
  PWeight ->
  Proposal (f Double)
hamiltonian x s n w = case checkHSettings s of
  Just err -> error err
  Nothing ->
    let desc = PDescription "Hamiltonian Monte Carlo (HMC)"
        dim = PSpecial (length x) 0.65
        ts = massesToTuningParameters (hMasses s)
        ps = hamiltonianSimple s
        p' = Proposal n desc dim w ps
        fT = defaultTuningFunction dim
        tS = hTune s
        fTs =
          if tS == HTuneMassesAndLeapfrog
            then computeAuxiliaryTuningParameters
            else \_ xs -> xs
     in case tS of
          HNoTune -> p' Nothing
          _ -> p' $ Just $ Tuner 1.0 fT ts fTs (hamiltonianSimpleWithTuningParameters s)
