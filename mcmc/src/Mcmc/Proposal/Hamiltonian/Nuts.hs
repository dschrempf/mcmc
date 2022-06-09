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
-- See Matthew D. Hoffman, Andrew Gelman (2014) The No-U-Turn Sampler:
-- Adaptively Setting Path Lengths in Hamiltonian Monte Carlo, Journal of
-- Machine Learning Research.
module Mcmc.Proposal.Hamiltonian.Nuts
  ( NTuningSpec,
    nTuningSpec,
    nuts,
  )
where

import Data.Bifunctor
import Mcmc.Proposal
import Mcmc.Proposal.Hamiltonian.Common
import Mcmc.Proposal.Hamiltonian.Internal
import Numeric.AD.Double
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import System.Random.MWC

-- Internal; Slice variable u.
type SliceVariable = Log Double

-- Internal; Forward is True.
type Direction = Bool

-- Internal; Doubling step number.
type DoublingStep = Int

-- Internal; The number of leapfrog steps within the slice ('n' in Algorithm 3).
type NStepsOk = Int

-- Internal; Well, that's fun, isn't it? Have a look at Algorithm 3.
-- cited above.
type BuildTreeReturnType = (Positions, Momenta, Positions, Momenta, Positions, NStepsOk)

-- Constant determining largest allowed leapfrog integration error. See
-- discussion around Equation (3).
deltaMax :: Log Double
deltaMax = Exp 1000

-- Second function in Algorithm 3.
buildTreeWith ::
  HData ->
  Target ->
  GenIO ->
  --
  Positions ->
  Momenta ->
  SliceVariable ->
  Direction ->
  DoublingStep ->
  LeapfrogScalingFactor ->
  IO (Maybe BuildTreeReturnType)
buildTreeWith hdata@(HData _ msInv) tfun g x p u v j e
  | j <= 0 =
      -- Move backwards or forwards?
      let e' = if v then e else negate e
       in case leapfrog tfun msInv 1 e' x p of
            Nothing -> pure Nothing
            Just (x', p', _, expEPot') ->
              if errorIsSmall
                then pure $ Just (x', p', x', p', x', n')
                else pure Nothing
              where
                expEKin' = exponentialKineticEnergy msInv p'
                n' = if u <= expEPot' * expEKin' then 1 else 0
                errorIsSmall = expEPot' * expEKin' > u / deltaMax
  -- Recursive case. This is complicated because the algorithm is written for an
  -- imperative language, and because we have two stacked monads.
  | otherwise = do
      mr <- buildTree x p u v (j - 1) e
      case mr of
        Nothing -> pure Nothing
        -- Here, 'm' stands for minus, and 'p' for plus.
        Just (xm, pm, xp, pp, x', n') -> do
          mr' <-
            if v
              then -- Forwards.
              do
                mr'' <- buildTree xp pp u v (j - 1) e
                case mr'' of
                  Nothing -> pure Nothing
                  Just (_, _, xp', pp', x'', n'') -> pure $ Just (xm, pm, xp', pp', x'', n'')
              else -- Backwards.
              do
                mr'' <- buildTree xm pm u v (j - 1) e
                case mr'' of
                  Nothing -> pure Nothing
                  Just (xm', pm', _, _, x'', n'') -> pure $ Just (xm', pm', xp, pp, x'', n'')
          case mr' of
            Nothing -> pure Nothing
            Just (xm'', pm'', xp'', pp'', x''', n''') -> do
              b <- uniform g :: IO Double
              let x'''' = if b < fromIntegral n''' / (fromIntegral $ n' + n''') then x''' else x'
                  n'''' = n' + n'''
                  -- Important: Check for U-turn. This formula differs from the
                  -- formula using indicator functions in Algorithm 3. However,
                  -- check Equation (4).
                  isUTurn = let dx = (xp'' - xm'') in (dx * pm'' < 0) || (dx * pp'' < 0)
              if isUTurn then pure Nothing else pure $ Just (xm'', pm'', xp'', pp'', x'''', n'''')
  where
    buildTree = buildTreeWith hdata tfun g

-- TODO (high): Tuning configuration.

-- TODO (high): Add tuning configuration to NTuningSpec.

-- | Complete tuning specification of the NUTS proposal.
--
-- Includes tuning parameters and tuning configuration.
data NTuningSpec = NTuningSpec
  { nMasses :: Masses,
    nLeapfrogScalingFactor :: LeapfrogScalingFactor
  }
  deriving (Show)

-- | See 'NTuningSpec'.
--
-- Return 'Left' if an error is found.
nTuningSpec :: Masses -> LeapfrogScalingFactor -> Either String NTuningSpec
nTuningSpec ms e
  | any (<= 0) diagonalMasses = eWith "Some diagonal entries of the mass matrix are zero or negative."
  | nrows /= ncols = eWith "Mass matrix is not square."
  | e <= 0 = eWith "Leapfrog scaling factor is zero or negative."
  | otherwise = Right $ NTuningSpec ms e
  where
    eWith m = Left $ "nTuningSpec: " <> m
    ms' = L.unSym ms
    diagonalMasses = L.toList $ L.takeDiag ms'
    nrows = L.rows ms'
    ncols = L.cols ms'

-- First function in Algorithm 3.
nutsPFunctionWithMemoizedCovariance ::
  NTuningSpec ->
  HSpec s ->
  HData ->
  (s Double -> Target) ->
  PFunction (s Double)
nutsPFunctionWithMemoizedCovariance nspec hspec hdata targetWith xComplete g = do
  p <- generateMomenta mu ms g
  uZeroOne <- uniform g :: IO Double
  -- NOTE (runtime): Here we need the target function value from the previous
  -- step. For now, I just recalculate the value, but this is, of course, slow!
  -- However, if other proposals have changed the state inbetween, we do need to
  -- recalculate this value.
  let x = toVec xComplete
      expEPot = fst $ target x
      expEKin = exponentialKineticEnergy msInv p
      uZeroOneL = Exp $ log uZeroOne
      u = expEPot * expEKin * uZeroOneL
  let -- Recursive case. This is complicated because the algorithm is written for an
      -- imperative language, and because we have two stacked monads.
      --
      -- Here, 'm' stands for minus, and 'p' for plus.
      go xm pm xp pp j y n isNewState = do
        v <- uniform g :: IO Direction
        mr' <-
          if v
            then -- Forwards.
            do
              mr <- buildTreeWith hdata target g xp pp u v j e
              case mr of
                Nothing -> pure Nothing
                Just (_, _, xp', pp', y', n') -> pure $ Just (xm, pm, xp', pp', y', n')
            else -- Backwards.
            do
              mr <- buildTreeWith hdata target g xm pm u v j e
              case mr of
                Nothing -> pure Nothing
                Just (xm', pm', _, _, y', n') -> pure $ Just (xm', pm', xp, pp, y', n')
        case mr' of
          Nothing -> pure (y, isNewState)
          Just (xm'', pm'', xp'', pp'', y'', n'') -> do
            let r = fromIntegral n'' / fromIntegral n
            isAccept <-
              if r > 1.0
                then pure True
                else do
                  b <- uniform g :: IO Double
                  pure $ b < r
            let (y''', isNewState') = if isAccept then (y'', True) else (y, isNewState)
                isUTurn = let dx = (xp'' - xm'') in (dx * pm'' < 0) || (dx * pp'' < 0)
            if isUTurn
              then pure (y''', isNewState')
              else go xm'' pm'' xp'' pp'' (j + 1) y''' (n + n'') isNewState'
  (x', isNew) <- go x p x p 0 x 1 False
  let r = if isNew then ForceAccept $ fromVec x' else ForceReject
  -- TODO @Dominik (high, feature): Acceptance counts.
  pure (r, Nothing)
  where
    (NTuningSpec ms e) = nspec
    (HSpec _ toVec fromVecWith) = hspec
    fromVec = fromVecWith xComplete
    (HData mu msInv) = hdata
    target = targetWith xComplete

nutsPFunction :: NTuningSpec -> HSpec s -> (s Double -> Target) -> PFunction (s Double)
nutsPFunction nspec hspec = nutsPFunctionWithMemoizedCovariance nspec hspec (getHData $ nMasses nspec)

-- | No U-turn Hamiltonian Monte Carlo sampler (NUTS).
nuts ::
  (Eq (s Double), Traversable s) =>
  NTuningSpec ->
  HSpec s ->
  HTarget s ->
  PName ->
  PWeight ->
  Proposal (s Double)
nuts nspec hspec htarget n w = case checkHSpecWith (nMasses nspec) hspec of
  Just err -> error err
  Nothing ->
    let -- Misc.
        desc = PDescription "No U-turn sampler (NUTS)"
        (HSpec sample toVec fromVec) = hspec
        dim = L.size $ toVec sample
        -- TODO (high): Desired acceptance ratio for NUTS?
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
        ps = nutsPFunction nspec hspec targetWith
        nutsWith = Proposal n desc PSlow pDim w ps
     in -- TODO (high): NUTS with tuning.
        nutsWith Nothing
