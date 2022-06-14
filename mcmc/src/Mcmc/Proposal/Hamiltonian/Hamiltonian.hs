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
-- The HMC proposal acts on 'Positions', a vector of floating point values. The
-- manipulated values can represent the complete state, or a subset of the
-- complete state. Functions converting the state to and from this vector have
-- to be provided; see 'HStructure'.
--
-- Even though the proposal may only act on a subset of the complete state, the
-- prior, likelihood, and Jacobian functions of the complete state have to be
-- provided; see 'HTarget'. This is because parameters not manipulated by the
-- HMC proposal still influence the prior, likelihood and Jacobian functions.
--
-- The points given above have implications on how the HMC proposal is handled:
-- Do not use 'liftProposalWith', 'liftProposal', or '(@~)' with the HMC
-- proposal; instead provide proper conversion functions with 'HStructure'.
--
-- The gradient of the log target function is calculated using automatic
-- differentiation; see the excellent
-- [ad](https://hackage.haskell.org/package/ad) package.
--
-- The desired acceptance rate is 0.65, although the dimension of the proposal
-- is high.
--
-- The speed of this proposal changes drastically with the leapfrog trajectory
-- length and the leapfrog scaling factor. Hence, the speed will change during
-- burn in.
--
-- References:
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
module Mcmc.Proposal.Hamiltonian.Hamiltonian
  ( -- * Hamiltonian Monte Carlo proposal
    HParams (..),
    defaultHParams,
    hamiltonian,
  )
where

import Data.Bifunctor
import Mcmc.Acceptance
import Mcmc.Algorithm.MHG
import Mcmc.Proposal
import Mcmc.Proposal.Hamiltonian.Common
import Mcmc.Proposal.Hamiltonian.Internal
import Mcmc.Proposal.Hamiltonian.Masses
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
-- - Estimate a reasonable leapfrog scaling factor using Algorithm 4 [4]. If all
--   fails, use 0.1.
--
-- - Leapfrog simulation length is set to 0.5.
--
-- - The mass matrix is set to the identity matrix.
defaultHParams :: HParams
defaultHParams = HParams Nothing Nothing Nothing

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

-- TODO @Dominik (high, feature): The expected acceptance counts should not be
-- calculated after burn in. Rather, the actual acceptance counts should be
-- reported. For this to work, the proposal needs to know if it is in "burn in
-- phase" or not.

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianPFunction'.
hamiltonianPFunction ::
  HParamsI ->
  HStructure s ->
  (s Double -> Target) ->
  PFunction (s Double)
hamiltonianPFunction hparamsi hstruct targetWith x g = do
  p <- generateMomenta mus ms g
  eRan <- uniformR (eL, eR) g
  -- NOTE: The NUTS paper does not sample l since l varies naturally because
  -- of epsilon. I still think it should vary because otherwise, there may be
  -- dragons due to periodicity.
  let lM = la / eRan
      lL = maximum [1 :: Int, floor $ 0.9 * lM]
      lR = maximum [lL, ceiling $ 1.1 * lM]
  lRan <- uniformR (lL, lR) g
  case leapfrog (targetWith x) msI lRan eRan q p of
    Nothing -> pure (ForceReject, Just $ AcceptanceCounts 0 100)
    -- Check if next state is accepted here, because the Jacobian is included in
    -- the target function. If not: pure (x, 0.0, 1.0).
    Just (q', p', prQ, prQ') -> do
      let -- Prior of momenta.
          prP = exponentialKineticEnergy msI p
          prP' = exponentialKineticEnergy msI p'
          r = prQ' * prP' / (prQ * prP)
      accept <- mhgAccept r g
      -- NOTE: For example, Neal page 12: In order for the Hamiltonian proposal
      -- to be in detailed balance, the momenta have to be negated before
      -- proposing the new value. That is, the negated momenta would guide the
      -- chain back to the previous state. However, we are only interested in
      -- the positions, and are not even storing the momenta.
      let pr = if accept then ForceAccept (fromVec x q') else ForceReject
          ar = exp $ ln r
          getCounts s = max 0 $ min 100 $ round $ s * 100
          ac =
            if ar >= 0
              then let cs = getCounts ar in AcceptanceCounts cs (100 - cs)
              else error $ "hamiltonianPFunction: Acceptance rate negative."
      pure (pr, Just ac)
  where
    (HParamsI e la ms _ _ msI mus) = hparamsi
    (HStructure _ toVec fromVec) = hstruct
    q = toVec x
    eL = 0.9 * e
    eR = 1.1 * e

-- | Hamiltonian Monte Carlo proposal.
--
-- The structure of the state is denoted as @s@.
--
-- May call 'error' during initialization.
hamiltonian ::
  Traversable s =>
  HParams ->
  HTuningConf ->
  HStructure s ->
  HTarget s ->
  PName ->
  PWeight ->
  Proposal (s Double)
hamiltonian hparams htconf hstruct htarget n w =
  let -- Misc.
      desc = PDescription "Hamiltonian Monte Carlo (HMC)"
      (HStructure sample toVec fromVec) = hstruct
      dim = L.size $ toVec sample
      -- See bottom of page 1615 in [4].
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
      (HParams mEps mLa mMs) = hparams
      hParamsI =
        either error id $
          hParamsIWith (targetWith sample) (toVec sample) mEps mLa mMs
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
