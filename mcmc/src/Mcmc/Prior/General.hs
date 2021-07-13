-- |
-- Module      :  Mcmc.Prior.General
-- Description :  Generalized prior functions for automatic differentiation
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Tue Jul 13 12:26:16 2021.
--
-- Generalized log prior functions. For the specialized versions, see
-- "Mcmc.Prior".
--
-- The generalized log prior functions act on 'RealFloat's to allow automatic
-- differentiation. Similar to the functions in "Mcmc.Prior", the results are in
-- __log domain__ although this is __not denoted by the type signatures__.
--
-- The generalized versions are tagged with a suffix @G@, and are expected to be
-- slower than the specialized ones in "Mcmc.Prior".
module Mcmc.Prior.General
  ( PriorFunctionG,

    -- * Improper priors
    noPriorG,
    greaterThanG,
    positiveG,
    lessThanG,
    negativeG,

    -- * Continuous priors
    exponentialG,
    gammaG,
    normalG,
    uniformG,
  )
where

import Mcmc.Internal.Gamma

-- | Generalized prior function.
type PriorFunctionG a = a -> a

-- | See 'Mcmc.Prior.noPrior'.
noPriorG :: RealFloat a => PriorFunctionG a
noPriorG = const 0.0
{-# SPECIALIZE noPriorG :: PriorFunctionG Double #-}

-- | See 'Mcmc.Prior.greaterThan'.
greaterThanG :: RealFloat a => a -> PriorFunctionG a
greaterThanG a x
  | x > a = 0
  | otherwise = - 1 / 0
{-# SPECIALIZE greaterThanG :: Double -> PriorFunctionG Double #-}

-- | See 'Mcmc.Prior.positive'.
positiveG :: RealFloat a => PriorFunctionG a
positiveG = greaterThanG 0
{-# SPECIALIZE positiveG :: PriorFunctionG Double #-}

-- | See 'Mcmc.Prior.lessThan'.
lessThanG :: RealFloat a => a -> PriorFunctionG a
lessThanG a x
  | x < a = 0.0
  | otherwise = - 1 / 0
{-# SPECIALIZE lessThanG :: Double -> PriorFunctionG Double #-}

-- | See 'Mcmc.Prior.negative'.
negativeG :: RealFloat a => PriorFunctionG a
negativeG = lessThanG 0
{-# SPECIALIZE negativeG :: PriorFunctionG Double #-}

-- | See 'Mcmc.Prior.exponential'.
exponentialG :: RealFloat a => a -> PriorFunctionG a
exponentialG l = (+ ll) . negate . (* l)
  where
    ll = log l
{-# SPECIALIZE exponentialG :: Double -> PriorFunctionG Double #-}

-- | See 'Mcmc.Prior.gamma'.
gammaG :: RealFloat a => a -> a -> PriorFunctionG a
gammaG k t x
  | x <= 0 = -1 / 0
  | otherwise = log x * (k - 1) - (x / t) - logGammaG k - log t * k
{-# SPECIALIZE gammaG :: Double -> Double -> PriorFunctionG Double #-}

mLnSqrt2Pi :: RealFloat a => a
mLnSqrt2Pi = 0.9189385332046727417803297364056176398613974736377834128171
{-# INLINE mLnSqrt2Pi #-}

-- | See 'Mcmc.Prior.normal'.
normalG :: RealFloat a => a -> a -> PriorFunctionG a
normalG m s x = (- xm * xm / (2 * s * s)) - denom
  where
    xm = x - m
    denom = mLnSqrt2Pi + log s
{-# SPECIALIZE normalG :: Double -> Double -> PriorFunctionG Double #-}

-- | See 'Mcmc.Prior.uniform'.
uniformG :: RealFloat a => a -> a -> PriorFunctionG a
uniformG a b x
  | x < a = - 1 / 0
  | x > b = - 1 / 0
  | otherwise = 0
