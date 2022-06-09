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
--   to be provided; see 'HSpec'.
--
-- - Even though the proposal may only act on a subset of the complete state,
--   the prior, likelihood, and Jacobian functions of the complete state have to
--   be provided; see 'HTarget'. This is because parameters not manipulated by
--   the HMC proposal still influence the prior, likelihood and Jacobian
--   functions.
--
-- - The points given above have implications on how the HMC proposal is
--   handled: Do not use 'liftProposalWith', 'liftProposal', or '(@~)' with the
--   HMC proposal; instead use the conversion functions in 'HSpec'.
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
  ( HTuningSpec,
    hTuningSpec,
    hamiltonian,
  )
where

import Data.Bifunctor
import qualified Data.Vector.Unboxed as VU
import Mcmc.Algorithm.MHG
import Mcmc.Proposal
import Mcmc.Proposal.Hamiltonian.Common
import Numeric.AD.Double
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import System.Random.MWC

-- | Complete tuning specification of the Hamilton Monte Carlo proposal.
--
-- Includes tuning parameters and tuning configuration.
data HTuningSpec = HTuningSpec
  { hMasses :: Masses,
    hLeapfrogTrajectoryLength :: LeapfrogTrajectoryLength,
    hLeapfrogScalingFactor :: LeapfrogScalingFactor,
    hTuningConf :: HTuningConf
  }
  deriving (Show)

-- | See 'HTuningSpec'.
--
-- Return 'Left' if an error is found.
hTuningSpec ::
  Masses ->
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  HTuningConf ->
  Either String HTuningSpec
hTuningSpec masses l eps c
  | any (<= 0) diagonalMasses = eWith "Some diagonal entries of the mass matrix are zero or negative."
  | nrows /= ncols = eWith "Mass matrix is not square."
  | l < 1 = eWith "Leapfrog trajectory length is zero or negative."
  | eps <= 0 = eWith "Leapfrog scaling factor is zero or negative."
  | otherwise = Right $ HTuningSpec masses l eps c
  where
    eWith m = Left $ "hTuningSpec: " <> m
    ms = L.unSym masses
    diagonalMasses = L.toList $ L.takeDiag ms
    nrows = L.rows ms
    ncols = L.cols ms

hTuningParametersToTuningSpec ::
  HTuningSpec ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String HTuningSpec
hTuningParametersToTuningSpec (HTuningSpec ms l e c) t ts
  | not nTsOK = Left "hTuningParametersToSettings: Auxiliary variables dimension mismatch."
  | otherwise = Right $ HTuningSpec msTuned lTuned eTuned c
  where
    d = L.rows $ L.unSym ms
    (HTuningConf tlf tms) = c
    nTsOK =
      let nTs = VU.length ts
       in case tms of
            HNoTuneMasses -> nTs == 0
            _ -> nTs == d * d
    msTuned = case tms of
      HNoTuneMasses -> ms
      _ -> tuningParametersToMasses d ts
    -- The larger epsilon, the larger the proposal step size and the lower the
    -- expected acceptance ratio.
    --
    -- Further, we roughly keep \( L * \epsilon = 1.0 \). The equation is not
    -- correct, because we pull L closer to the original value to keep the
    -- runtime somewhat acceptable.
    (lTuned, eTuned) = case tlf of
      HNoTuneLeapfrog -> (l, e)
      HTuneLeapfrog -> (ceiling $ fromIntegral l / (t ** 0.8) :: Int, t * e)

hamiltonianProposeWithTuningParameters ::
  Traversable s =>
  HTuningSpec ->
  HSpec s ->
  (s Double -> Target) ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (Propose (s Double))
hamiltonianProposeWithTuningParameters tspec hspec targetWith t ts = do
  tspec' <- hTuningParametersToTuningSpec tspec t ts
  pure $ hamiltonianPropose tspec' hspec targetWith

-- The inverted covariance matrix and the log determinant of the covariance
-- matrix are calculated by 'hamiltonianPropose'.
hamiltonianProposeWithMemoizedCovariance ::
  Traversable s =>
  HTuningSpec ->
  HSpec s ->
  HData ->
  (s Double -> Target) ->
  Propose (s Double)
hamiltonianProposeWithMemoizedCovariance tspec hspec dt targetWith x g = do
  phi <- generateMomenta mu masses g
  lRan <- uniformR (lL, lR) g
  eRan <- uniformR (eL, eR) g
  case leapfrog (targetWith x) massesInv lRan eRan theta phi of
    -- TODO @Dominik (high, feature): Acceptance counts.
    Nothing -> pure (ForceReject, Nothing)
    -- Check if next state is accepted here, because the Jacobian is included in
    -- the target function. If not: pure (x, 0.0, 1.0).
    Just (theta', phi', prTheta, prTheta') -> do
      let -- Prior of momenta.
          prPhi = exponentialKineticEnergy massesInv phi
          prPhi' = exponentialKineticEnergy massesInv phi'
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
    (HTuningSpec masses l e _) = tspec
    (HSpec _ toVec fromVec) = hspec
    theta = toVec x
    lL = maximum [1 :: Int, floor $ (0.8 :: Double) * fromIntegral l]
    lR = maximum [lL, ceiling $ (1.2 :: Double) * fromIntegral l]
    eL = 0.8 * e
    eR = 1.2 * e
    (HData mu massesInv) = dt

hamiltonianPropose ::
  Traversable s =>
  HTuningSpec ->
  HSpec s ->
  (s Double -> Target) ->
  Propose (s Double)
hamiltonianPropose tspec hspec = hamiltonianProposeWithMemoizedCovariance tspec hspec (getHData $ hMasses tspec)

hGetTuningFunction :: Int -> (a -> Positions) -> HTuningConf -> Maybe (TuningFunction a)
hGetTuningFunction n toVec (HTuningConf l m) = case (l, m) of
  (HNoTuneLeapfrog, HNoTuneMasses) -> Nothing
  (HNoTuneLeapfrog, HTuneDiagonalMassesOnly) -> Just $ tuningFunctionOnlyAux td
  (HNoTuneLeapfrog, HTuneAllMasses) -> Just $ tuningFunctionOnlyAux ta
  (HTuneLeapfrog, HNoTuneMasses) -> Just tuningFunction
  (HTuneLeapfrog, HTuneDiagonalMassesOnly) -> Just $ tuningFunctionWithAux td
  (HTuneLeapfrog, HTuneAllMasses) -> Just $ tuningFunctionWithAux ta
  where
    td = tuneDiagonalMassesOnly n toVec
    ta = tuneAllMasses n toVec

-- | Hamiltonian Monte Carlo proposal.
hamiltonian ::
  (Eq (s Double), Traversable s) =>
  HTuningSpec ->
  HSpec s ->
  HTarget s ->
  PName ->
  PWeight ->
  Proposal (s Double)
hamiltonian tspec hspec htarget n w = case checkHSpecWith (hMasses tspec) hspec of
  Just err -> error err
  Nothing ->
    let -- Misc.
        desc = PDescription "Hamiltonian Monte Carlo (HMC)"
        (HSpec sample toVec fromVec) = hspec
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
        ps = hamiltonianPropose tspec hspec targetWith
        hamiltonianWith = Proposal n desc PSlow pDim w ps
        -- Tuning.
        tconf@(HTuningConf _ tms) = hTuningConf tspec
        ts = case tms of
          HNoTuneMasses -> VU.empty
          _ -> massesToTuningParameters $ hMasses tspec
        tuner = do
          tfun <- hGetTuningFunction dim toVec tconf
          pure $ Tuner 1.0 ts tfun (hamiltonianProposeWithTuningParameters tspec hspec targetWith)
     in hamiltonianWith tuner
