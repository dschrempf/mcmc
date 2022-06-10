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
-- - [4] Matthew D. Hoffman, Andrew Gelman (2014) The No-U-Turn Sampler:
--   Adaptively Setting Path Lengths in Hamiltonian Monte Carlo, Journal of
--   Machine Learning Research.
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
import Mcmc.Acceptance
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
-- If a parameter is 'Nothing', a default value is used (see 'defaultHParams').
data HParams = HParams
  { hLeapfrogScalingFactor :: Maybe LeapfrogScalingFactor,
    hLeapfrogSimulationLength :: Maybe LeapfrogSimulationLength,
    hMasses :: Maybe Masses
  }
  deriving (Show)

-- | Default parameters.
--
-- - Estimate a reasonable leapfrog scaling factor using Algorithm 4 in Matthew
--   D. Hoffman, Andrew Gelman (2014) The No-U-Turn Sampler: Adaptively Setting
--   Path Lengths in Hamiltonian Monte Carlo, Journal of Machine Learning
--   Research.
--
-- - Leapfrog simulation length is set to 0.5.
--
-- - The mass matrix is set to the identity matrix.
defaultHParams :: HParams
defaultHParams = HParams Nothing Nothing Nothing

-- Check 'HParamsI; set parameters; compute 'HData'.
--
-- NOTE: This function may sample momenta to find a reasonable leapfrog scaling
-- factor. Hence, it should be run in the IO monad. However, I opt for a pure
-- setting here and use a fixed seed.
fromHParams :: Target -> Positions -> HParams -> Either String HParamsI
fromHParams htarget p (HParams mEps mLa mMs) = do
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
  la <- case mLa of
    Nothing -> Right 0.5
    Just l
      | l <= 0 -> eWith "Leapfrog simulation length is zero or negative."
      | otherwise -> Right l
  eps <- case mEps of
    Nothing -> Right $ runST $ do
      g <- create
      findReasonableEpsilon htarget ms p g
    Just e
      | e <= 0 -> eWith "Leapfrog scaling factor is zero or negative."
      | otherwise -> Right e
  pure $ hParamsIWith eps la ms
  where
    eWith m = Left $ "fromHParams: " <> m

hamiltonianPFunctionWithTuningParameters ::
  Traversable s =>
  Dimension ->
  HStructure s ->
  (s Double -> Target) ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (PFunction (s Double))
hamiltonianPFunctionWithTuningParameters d hstruct targetWith _ ts = do
  hParamsI <- fromAuxiliaryTuningParameters d ts
  pure $ hamiltonianPFunction hParamsI hstruct targetWith

-- TODO @Dominik (high, issue): Acceptance counts. How to combine with values
-- reported here and from the NUTS sampler.

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianPFunction'.
hamiltonianPFunctionWithMemoizedCovariance ::
  Traversable s =>
  HParamsI ->
  HStructure s ->
  (s Double -> Target) ->
  PFunction (s Double)
hamiltonianPFunctionWithMemoizedCovariance tspec hstruct targetWith x g = do
  phi <- generateMomenta mu ms g
  eRan <- uniformR (eL, eR) g
  -- NOTE: The NUTS paper does not sample l since l varies naturally because
  -- of epsilon. I still think it should vary because otherwise, there may be
  -- dragons due to periodicity.
  let lM = la / eRan
      lL = maximum [1 :: Int, floor $ 0.9 * lM]
      lR = maximum [lL, ceiling $ 1.1 * lM]
  lRan <- uniformR (lL, lR) g
  case leapfrog (targetWith x) msInv lRan eRan theta phi of
    Nothing -> pure (ForceReject, Just $ AcceptanceCounts 0 100)
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
          ar = exp $ ln r
          ac =
            if ar >= 0
              then let a = max 100 (round (ar * 100)) in AcceptanceCounts a (100 - a)
              else error $ "hamiltonianPFunctionWithMemoizedCovariance: Acceptance rate negative." <> show ar
      pure (pr, Just ac)
  where
    (HParamsI e la ms _ _ hdata) = tspec
    (HData mu msInv) = hdata
    -- TODO: The sample should not be in HStructure.
    (HStructure _ toVec fromVec) = hstruct
    theta = toVec x
    eL = 0.9 * e
    eR = 1.1 * e

hamiltonianPFunction ::
  Traversable s =>
  HParamsI ->
  HStructure s ->
  (s Double -> Target) ->
  PFunction (s Double)
hamiltonianPFunction hparams hstruct = hamiltonianPFunctionWithMemoizedCovariance hparams hstruct

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
      -- See bottom of page 1615 in Matthew D. Hoffman, Andrew Gelman (2014) The
      -- No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte
      -- Carlo, Journal of Machine Learning Research.
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
      hParamsI = either error id $ fromHParams (targetWith sample) (toVec sample) hparams
      ps = hamiltonianPFunction hParamsI hstruct targetWith
      hamiltonianWith = Proposal n desc PSlow pDim w ps
      -- Tuning.
      ts = toAuxiliaryTuningParameters hParamsI
      tuner = do
        tfun <- hTuningFunctionWith dim toVec htconf
        let pfun = hamiltonianPFunctionWithTuningParameters dim hstruct targetWith
        pure $ Tuner 1.0 ts tfun pfun
   in case checkHStructureWith (hpsMasses hParamsI) hstruct of
        Just err -> error err
        Nothing -> hamiltonianWith tuner
