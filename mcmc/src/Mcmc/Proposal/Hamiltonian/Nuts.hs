-- |
-- Module      :  Mcmc.Proposal.Hamiltonian.Nuts
-- Description :  No-U-Turn sampler (NUTS)
-- Copyright   :  2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri May 27 09:58:23 2022.
--
-- For a general introduction to Hamiltonian proposals, see
-- "Mcmc.Proposal.Hamiltonian.Hamiltonian".
--
-- This module implements the No-U-Turn Sampler (NUTS), as described in [4].
--
-- Work in progress.
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
module Mcmc.Proposal.Hamiltonian.Nuts
  ( NParams (..),
    defaultNParams,
    nuts,
  )
where

import Data.Bifunctor
import Mcmc.Acceptance
import Mcmc.Proposal
import Mcmc.Proposal.Hamiltonian.Common
import Mcmc.Proposal.Hamiltonian.Internal
import Mcmc.Proposal.Hamiltonian.Masses
import Numeric.AD.Double
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import System.Random.Stateful

-- Internal; Slice variable 'u'.
type SliceVariable = Log Double

-- Internal; Forward is True.
type Direction = Bool

-- Internal; Doubling step number 'j'.
type DoublingStep = Int

-- Internal; Number of leapfrog steps within the slice 'n'.
type NStepsOk = Int

-- Internal; Estimated acceptance rate \(\alpha\)'.
type Alpha = Log Double

-- Internal; Number of accepted steps.
type NAlpha = Int

-- Internal; Well, that's fun, isn't it? Have a look at Algorithm 3 in [4].
type BuildTreeReturnType = (Positions, Momenta, Positions, Momenta, Positions, NStepsOk, Alpha, NAlpha)

-- Constant determining largest allowed leapfrog integration error. See
-- discussion around Equation (3) in [4].
deltaMax :: Log Double
deltaMax = Exp 1000

-- Second function in Algorithm 3 and Algorithm 6, respectively in [4].
buildTreeWith ::
  -- The exponent of the total energy of the starting state is used to
  -- calculate the expected acceptance rate 'Alpha'.
  Log Double ->
  MassesI ->
  Target ->
  IOGenM StdGen ->
  --
  Positions ->
  Momenta ->
  SliceVariable ->
  Direction ->
  DoublingStep ->
  LeapfrogScalingFactor ->
  IO (Maybe BuildTreeReturnType)
buildTreeWith expETot0 msI tfun g q p u v j e
  | j <= 0 =
      -- Move backwards or forwards?
      let e' = if v then e else negate e
       in case leapfrog tfun msI 1 e' q p of
            Nothing -> pure Nothing
            Just (q', p', _, expEPot') ->
              if errorIsSmall
                then pure $ Just (q', p', q', p', q', n', alpha, 1)
                else pure Nothing
              where
                expEKin' = exponentialKineticEnergy msI p'
                expETot' = expEPot' * expEKin'
                n' = if u <= expEPot' * expEKin' then 1 else 0
                errorIsSmall = u < deltaMax * expETot'
                alpha' = expETot' / expETot0
                -- Limit expected acceptance rate between 0 and 1.
                alpha = max 0 $ min 1 alpha'
  -- Recursive case. This is complicated because the algorithm is written for an
  -- imperative language, and because we have two stacked monads.
  | otherwise = do
      mr <- buildTree q p u v (j - 1) e
      case mr of
        Nothing -> pure Nothing
        -- Here, the suffixes 'm' and 'p' stand for minus and plus, respectively.
        Just (qm, pm, qp, pp, q', n', a', na') -> do
          mr' <-
            if v
              then -- Forwards.
              do
                mr'' <- buildTree qp pp u v (j - 1) e
                case mr'' of
                  Nothing -> pure Nothing
                  Just (_, _, qp', pp', q'', n'', a'', na'') ->
                    pure $ Just (qm, pm, qp', pp', q'', n'', a'', na'')
              else -- Backwards.
              do
                mr'' <- buildTree qm pm u v (j - 1) e
                case mr'' of
                  Nothing -> pure Nothing
                  Just (qm', pm', _, _, q'', n'', a'', na'') ->
                    pure $ Just (qm', pm', qp, pp, q'', n'', a'', na'')
          case mr' of
            Nothing -> pure Nothing
            Just (qm'', pm'', qp'', pp'', q''', n''', a''', na''') -> do
              b <- uniformRM (0, 1) g :: IO Double
              let n'''' = n' + n'''
                  q'''' = if b < fromIntegral n''' / fromIntegral n'''' then q''' else q'
                  a'''' = a' + a'''
                  na'''' = na' + na'''
                  -- Important: Check for U-turn. This formula differs from the
                  -- formula using indicator functions in Algorithm 3. However,
                  -- check Equation (4).
                  isUTurn = let dq = (qp'' - qm'') in (dq * pm'' < 0) || (dq * pp'' < 0)
              if isUTurn
                then pure Nothing
                else pure $ Just (qm'', pm'', qp'', pp'', q'''', n'''', a'''', na'''')
  where
    buildTree = buildTreeWith expETot0 msI tfun g

-- | Parameters of the NUTS proposal.
--
-- Includes tuning parameters and tuning configuration.
data NParams = NParams
  { nLeapfrogScalingFactor :: Maybe LeapfrogScalingFactor,
    nMasses :: Maybe Masses
  }
  deriving (Show)

-- | Default parameters.
--
-- - Estimate a reasonable leapfrog scaling factor using Algorithm 4 [4]. If all
--   fails, use 0.1.
--
-- - The mass matrix is set to the identity matrix.
defaultNParams :: NParams
defaultNParams = NParams Nothing Nothing

nutsPFunctionWithTuningParameters ::
  (Traversable s) =>
  Dimension ->
  HStructure s ->
  (s Double -> Target) ->
  TuningParameter ->
  AuxiliaryTuningParameters ->
  Either String (PFunction (s Double))
nutsPFunctionWithTuningParameters d hstruct targetWith _ ts = do
  hParamsI <- fromAuxiliaryTuningParameters d ts
  pure $ nutsPFunction hParamsI hstruct targetWith

data IsNew
  = Old
  | OldWith {_acceptanceRatesOld :: AcceptanceRates}
  | NewWith {_acceptanceRatesNew :: AcceptanceRates}

-- First function in Algorithm 3.
nutsPFunction ::
  HParamsI ->
  HStructure s ->
  (s Double -> Target) ->
  PFunction (s Double)
nutsPFunction hparamsi hstruct targetWith x g = do
  p <- generateMomenta mus ms g
  uZeroOne <- uniformRM (0, 1) g :: IO Double
  -- NOTE (runtime): Here we need the target function value from the previous
  -- step. For now, I just recalculate the value, but this is, of course, slow!
  -- However, if other proposals have changed the state inbetween, we do need to
  -- recalculate this value.
  let q = toVec x
      expEPot = fst $ target q
      expEKin = exponentialKineticEnergy msI p
      expETot = expEPot * expEKin
      uZeroOneL = Exp $ log uZeroOne
      u = expETot * uZeroOneL
  let -- Recursive case. This is complicated because the algorithm is written for an
      -- imperative language, and because we have two stacked monads.
      --
      -- Here, the suffixes 'm' and 'p' stand for minus and plus, respectively.
      go qm pm qp pp j y n isNew = do
        v <- uniformM g :: IO Direction
        mr' <-
          if v
            then -- Forwards.
            do
              mr <- buildTreeWith expETot msI target g qp pp u v j e
              case mr of
                Nothing -> pure Nothing
                Just (_, _, qp', pp', y', n', a, na) -> pure $ Just (qm, pm, qp', pp', y', n', a, na)
            else -- Backwards.
            do
              mr <- buildTreeWith expETot msI target g qm pm u v j e
              case mr of
                Nothing -> pure Nothing
                Just (qm', pm', _, _, y', n', a, na) -> pure $ Just (qm', pm', qp, pp, y', n', a, na)
        case mr' of
          Nothing -> pure (y, isNew)
          Just (qm'', pm'', qp'', pp'', y'', n'', a, na) -> do
            let r = fromIntegral n'' / fromIntegral n :: Double
                -- Individual expected acceptance rates are limited in
                -- 'buildTreeWith'.
                ar = max 0 $ exp $ ln a
                ars = AcceptanceRates ar na
            isAccept <-
              if r > 1.0
                then pure True
                else do
                  b <- uniformRM (0, 1) g
                  pure $ b < r
            let (y''', isNew') = if isAccept then (y'', NewWith ars) else (y, OldWith ars)
                isUTurn = let dq = (qp'' - qm'') in (dq * pm'' < 0) || (dq * pp'' < 0)
            if isUTurn
              then pure (y''', isNew')
              else go qm'' pm'' qp'' pp'' (j + 1) y''' (n + n'') isNew'
  (x', isNew) <- go q p q p 0 q 1 Old
  pure $ case isNew of
    -- NOTE: I am not sure if it is correct to set the expected acceptance rate
    -- to 0 when the no u-turn leapfrog integrator fails.
    Old -> (ForceReject, Just $ AcceptanceRates 0 1)
    OldWith ac -> (ForceReject, Just ac)
    NewWith ac -> (ForceAccept $ fromVec x', Just ac)
  where
    (HParamsI e _ ms _ _ msI mus) = hparamsi
    (HStructure _ toVec fromVecWith) = hstruct
    fromVec = fromVecWith x
    target = targetWith x

-- | No U-turn Hamiltonian Monte Carlo sampler (NUTS).
--
-- The structure of the state is denoted as @s@.
--
-- May call 'error' during initialization.
nuts ::
  (Traversable s) =>
  NParams ->
  HTuningConf ->
  HStructure s ->
  HTarget s ->
  PName ->
  PWeight ->
  Proposal (s Double)
nuts nparams htconf hstruct htarget n w =
  let -- Misc.
      desc = PDescription "No U-turn sampler (NUTS)"
      (HStructure sample toVec fromVec) = hstruct
      dim = L.size $ toVec sample
      -- See bottom of page 1616 in [4].
      pDim = PSpecial dim 0.6
      -- Vectorize and derive the target function.
      (HTarget mPrF lhF mJcF) = htarget
      tF y = case (mPrF, mJcF) of
        (Nothing, Nothing) -> lhF y
        (Just prF, Nothing) -> prF y * lhF y
        (Nothing, Just jcF) -> lhF y * jcF y
        (Just prF, Just jcF) -> prF y * lhF y * jcF y
      tFnG = grad' (ln . tF)
      targetWith x = bimap Exp toVec . tFnG . fromVec x
      (NParams mEps mMs) = nparams
      hParamsI =
        either error id $
          hParamsIWith (targetWith sample) (toVec sample) mEps Nothing mMs
      ps = nutsPFunction hParamsI hstruct targetWith
      nutsWith = Proposal n desc PSlow pDim w ps
      -- Tuning.
      ts = toAuxiliaryTuningParameters hParamsI
      tuner = do
        tfun <- hTuningFunctionWith dim toVec htconf
        let pfun = nutsPFunctionWithTuningParameters dim hstruct targetWith
        pure $ Tuner 1.0 ts True True tfun pfun
   in case checkHStructureWith (hpsMasses hParamsI) hstruct of
        Just err -> error err
        Nothing -> nutsWith tuner
