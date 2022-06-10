{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      :  Mcmc.Proposal.Hamiltonian.Hamiltonian
-- Description :  Hamiltonian Monte Carlo proposal
-- Copyright   :  2021 Dominik Schrempf
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
-- Notes on implementation:
--
-- - The HMC proposal acts on 'Positions', a vector of floating point values.
--   The manipulated values can represent the complete state, or a subset of the
--   complete state. Functions converting the state to and from this vector have
--   to be provided; see 'HStructure'.
--
-- - Even though the proposal may only act on a subset of the complete state,
--   the prior, likelihood, and Jacobian functions of the complete state have to
--   be provided; see 'HTarget'. This is because parameters not manipulated by
--   the HMC proposal still influence the prior, likelihood and Jacobian
--   functions.
--
-- - The points given above have implications on how the HMC proposal is
--   handled: Do not use 'liftProposalWith', 'liftProposal', or '(@~)' with the
--   HMC proposal; instead use the conversion functions in 'HStructure'.
--
-- - The gradient of the log target function is calculated using automatic
--   differentiation.
--
-- - The desired acceptance rate is 0.65, although the dimension of the proposal
--   is high.
--
-- - The speed of this proposal changes drastically with the leapfrog trajectory
--   length and the leapfrog scaling factor. Hence, the speed will change during
--   burn in.
module Mcmc.Proposal.Hamiltonian.Hamiltonian
  ( -- * Hamiltonian Monte Carlo proposal
    HParams (..),
    defaultHParams,
    hamiltonian,
  )
where

import Control.Monad
import Control.Monad.ST
import Data.Bifunctor
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Mcmc.Algorithm.MHG
import Mcmc.Proposal
import Mcmc.Proposal.Hamiltonian.Common
import Mcmc.Proposal.Hamiltonian.Internal
import Numeric.AD.Double
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import System.Random.MWC

-- | Parameters of the Hamilton Monte Carlo proposal.
--
-- If a parameter is 'Nothing', the value in 'defaultHParams' will be used.
data HParams = HParams
  { hMasses :: Maybe Masses,
    hLeapfrogTrajectoryLength :: Maybe LeapfrogTrajectoryLength,
    hLeapfrogScalingFactor :: Maybe LeapfrogScalingFactor
  }
  deriving (Show)

-- | Default parameters.
--
-- - Mass matrix is identity matrix.
--
-- - Leapfrog trajectory length is 10.
--
-- - A resonable leapfrog scaling factor will be estimated.
defaultHParams :: HParams
defaultHParams = HParams Nothing Nothing Nothing

-- Internal. We need two data types here. 'HParams' is exposed to the
-- user. This one is exposed to the algorithm, and has all values set.
data HParamsSet = HParamsSet
  { hMasses' :: Masses,
    hLeapfrogTrajectoryLength' :: LeapfrogTrajectoryLength,
    hLeapfrogScalingFactor' :: LeapfrogScalingFactor,
    hData :: HData
  }
  deriving (Show)

-- Check 'HParams'; set parameters; compute 'HData'.
--
-- NOTE: This function may sample momenta to find a reasonable leapfrog scaling
-- factor. Hence, it should be run in the IO monad. However, I opt for a pure
-- setting here and use a fixed seed.
fromHParams :: Target -> Positions -> HParams -> Either String HParamsSet
fromHParams htarget p (HParams mMs mL mE) = do
  d <- case VS.length p of
    0 -> eWith "Empty position vector."
    d -> Right d
  let defMs = L.trustSym $ L.ident d
  ms <- case mMs of
    Nothing -> Right defMs
    Just ms -> do
      let ms' = L.unSym ms
          diagonalMs = L.toList $ L.takeDiag ms'
      when (any (<= 0) diagonalMs) $ eWith "Some diagonal masses are zero or negative."
      let nrows = L.rows ms'
          ncols = L.cols ms'
      when (nrows /= ncols) $ eWith "Mass matrix is not square."
      Right ms
  let hdata = getHData ms
  l <- case mL of
    Nothing -> Right 10
    Just l
      | l <= 0 -> eWith "Leapfrog trajectory length is zero or negative."
      | otherwise -> Right l
  e <- case mE of
    Nothing -> Right $ runST $ do
      g <- create
      findReasonableEpsilon htarget ms hdata p g
    Just e
      | e <= 0 -> eWith "Leapfrog scaling factor is zero or negative."
      | otherwise -> Right e
  pure $ HParamsSet ms l e hdata
  where
    eWith m = Left $ "fromHParams: " <> m

hTuningParametersToHParamsSet ::
  HParamsSet ->
  HTuningConf ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String HParamsSet
hTuningParametersToHParamsSet (HParamsSet ms l e hdata) c t ts
  | not nTsOK = Left "hTuningParametersToHParamsSet: Auxiliary variables dimension mismatch."
  | otherwise = Right $ HParamsSet msTuned lTuned eTuned hdataTuned
  where
    d = L.rows $ L.unSym ms
    (HTuningConf tlf tms) = c
    nTsOK =
      let nTs = VU.length ts
       in case tms of
            HNoTuneMasses -> nTs == 0
            _ -> nTs == d * d
    (msTuned, hdataTuned) = case tms of
      HNoTuneMasses -> (ms, hdata)
      _ -> let ms' = tuningParametersToMasses d ts in (ms', getHData ms')
    -- The larger epsilon, the larger the proposal step size and the lower the
    -- expected acceptance ratio.
    --
    -- Further, we roughly keep \( L * \epsilon = 1.0 \). The equation is not
    -- correct, because we pull L closer to the original value to keep the
    -- runtime somewhat acceptable.
    (lTuned, eTuned) = case tlf of
      HNoTuneLeapfrog -> (l, e)
      HTuneLeapfrog -> (ceiling $ fromIntegral l / (t ** 0.8) :: Int, t * e)

hamiltonianPFunctionWithTuningParameters ::
  Traversable s =>
  HParamsSet ->
  HTuningConf ->
  HStructure s ->
  (s Double -> Target) ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (PFunction (s Double))
hamiltonianPFunctionWithTuningParameters tspec htconf hstruct targetWith t ts = do
  tspec' <- hTuningParametersToHParamsSet tspec htconf t ts
  pure $ hamiltonianPFunction tspec' hstruct targetWith

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianPFunction'.
hamiltonianPFunctionWithMemoizedCovariance ::
  Traversable s =>
  HParamsSet ->
  HStructure s ->
  (s Double -> Target) ->
  PFunction (s Double)
hamiltonianPFunctionWithMemoizedCovariance tspec hstruct targetWith x g = do
  phi <- generateMomenta mu ms g
  lRan <- uniformR (lL, lR) g
  eRan <- uniformR (eL, eR) g
  case leapfrog (targetWith x) msInv lRan eRan theta phi of
    -- TODO @Dominik (high, feature): Acceptance counts.
    Nothing -> pure (ForceReject, Nothing)
    -- Check if next state is accepted here, because the Jacobian is included in
    -- the target function. If not: pure (x, 0.0, 1.0).
    Just (theta', phi', prTheta, prTheta') -> do
      let -- Prior of momenta.
          prPhi = exponentialKineticEnergy msInv phi
          prPhi' = exponentialKineticEnergy msInv phi'
          r = prTheta' * prPhi' / (prTheta * prPhi)
      accept <- mhgAccept r g
      -- NOTE: For example, Neal page 12: In order for the Hamiltonian proposal
      -- to be in detailed balance, the momenta have to be negated before
      -- proposing the new value. That is, the negated momenta would guide the
      -- chain back to the previous state. However, we are only interested in
      -- the positions, and are not even storing the momenta.
      let pr = if accept then ForceAccept (fromVec x theta') else ForceReject
      -- TODO @Dominik (high, feature): Acceptance counts.
      pure (pr, Nothing)
  where
    (HParamsSet ms l e hdata) = tspec
    (HData mu msInv) = hdata
    (HStructure _ toVec fromVec) = hstruct
    theta = toVec x
    lL = maximum [1 :: Int, floor $ (0.8 :: Double) * fromIntegral l]
    lR = maximum [lL, ceiling $ (1.2 :: Double) * fromIntegral l]
    eL = 0.8 * e
    eR = 1.2 * e

hamiltonianPFunction ::
  Traversable s =>
  HParamsSet ->
  HStructure s ->
  (s Double -> Target) ->
  PFunction (s Double)
hamiltonianPFunction hparams hstruct = hamiltonianPFunctionWithMemoizedCovariance hparams hstruct

hGetTuningFunction :: Int -> (a -> Positions) -> HTuningConf -> Maybe (TuningFunction a)
hGetTuningFunction n toVec (HTuningConf l m) = case (l, m) of
  (HNoTuneLeapfrog, HNoTuneMasses) -> Nothing
  (HNoTuneLeapfrog, HTuneDiagonalMassesOnly) -> Just $ tuningFunctionOnlyAux td
  (HNoTuneLeapfrog, HTuneAllMasses) -> Just $ tuningFunctionOnlyAux ta
  (HTuneLeapfrog, HNoTuneMasses) -> Just tuningFunction
  (HTuneLeapfrog, HTuneDiagonalMassesOnly) -> Just $ tuningFunctionWithAux td
  (HTuneLeapfrog, HTuneAllMasses) -> Just $ tuningFunctionWithAux ta
  where
    td = const (tuneDiagonalMassesOnly n toVec)
    ta = const (tuneAllMasses n toVec)

-- | Hamiltonian Monte Carlo proposal.
--
-- May call 'error' during initialization.
hamiltonian ::
  (Eq (s Double), Traversable s) =>
  HParams ->
  HTuningConf ->
  HStructure s ->
  HTarget s ->
  PName ->
  PWeight ->
  -- NOTE: Should probably be an Either.
  Proposal (s Double)
hamiltonian hparams htconf hstruct htarget n w =
  let -- Misc.
      desc = PDescription "Hamiltonian Monte Carlo (HMC)"
      (HStructure sample toVec fromVec) = hstruct
      dim = L.size $ toVec sample
      pDim = PSpecial dim 0.65
      -- Vectorize and derive the target function.
      (HTarget mPrF lhF mJcF) = htarget
      tF y = case (mPrF, mJcF) of
        (Nothing, Nothing) -> lhF y
        (Just prF, Nothing) -> prF y * lhF y
        (Nothing, Just jcF) -> lhF y * jcF y
        (Just prF, Just jcF) -> prF y * lhF y * jcF y
      tFnG = grad' (ln . tF)
      targetWith x = bimap Exp toVec . tFnG . fromVec x
      hparams' = either error id $ fromHParams (targetWith sample) (toVec sample) hparams
      ps = hamiltonianPFunction hparams' hstruct targetWith
      hamiltonianWith = Proposal n desc PSlow pDim w ps
      -- Tuning.
      (HTuningConf _ tms) = htconf
      ts = case tms of
        HNoTuneMasses -> VU.empty
        _ -> massesToTuningParameters $ hMasses' hparams'
      tuner = do
        tfun <- hGetTuningFunction dim toVec htconf
        pure $ Tuner 1.0 ts tfun (hamiltonianPFunctionWithTuningParameters hparams' htconf hstruct targetWith)
   in case checkHStructureWith (hMasses' hparams') hstruct of
        Just err -> error err
        Nothing -> hamiltonianWith tuner
