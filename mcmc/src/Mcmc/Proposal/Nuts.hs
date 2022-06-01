-- |
-- Module      :  Mcmc.Proposal.Nuts
-- Description :  No-U-Turn sampler (NUTS)
-- Copyright   :  (c) 2022 Dominik Schrempf
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
module Mcmc.Proposal.Nuts
  (
  )
where

import Mcmc.Proposal.Hamiltonian
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
deltaMax :: Double
deltaMax = 1000

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
buildTreeWith hdata@(HData mu msInv logDetMs) tfun g x p u v j e
  | j <= 0 =
      -- Move backwards or forwards?
      let e' = if v then e else negate e
       in case leapfrog tfun msInv 1 e' x p of
            Nothing -> pure Nothing
            Just (x', p', _, ePot') ->
              if errorIsLarge
                then pure Nothing
                else pure $ Just (x', p', x', p', x', n')
              where
                eKin' = logDensityMultivariateNormal mu msInv logDetMs p'
                n' = if u <= ePot' * eKin' then 1 else 0
                errorIsLarge = ln ePot' + ln eKin' > ln u - deltaMax
  | otherwise = do
      mr <- buildTreeWith hdata tfun g x p u v (j - 1) e
      case mr of
        Nothing -> pure Nothing
        -- Here, 'm' stands for minus, and 'p' for plus.
        Just (xm, pm, xp, pp, x', n') -> do
          mr' <-
            if v
              then -- Forwards.
              do
                mr'' <- buildTreeWith hdata tfun g xm pm u v (j - 1) e
                case mr'' of
                  Nothing -> pure Nothing
                  Just (_, _, xp', pp', x'', n'') -> pure $ Just (xm, pm, xp', pp', x'', n'')
              else -- Backwards.
              do
                mr'' <- buildTreeWith hdata tfun g xp pp u v (j - 1) e
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
