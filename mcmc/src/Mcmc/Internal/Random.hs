{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      :  Mcmc.Internal.Random
-- Description :  Tools for random calculations
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Nov 25 07:14:52 2020.
module Mcmc.Internal.Random
  ( splitGen,
    saveGen,
    loadGen,
  )
where

import Control.Monad
import Control.Monad.Primitive
import qualified Data.Vector.Unboxed as V
import Data.Word
import System.Random.MWC
import System.IO.Unsafe

-- | Split a generator.
--
-- Splitting an MWC pseudo number generator is not good practice. However, I
-- have to go with this solution for now, and wait for proper support of
-- spittable pseudo random number generators such as @splitmix@.
splitGen :: PrimMonad m => Int -> Gen (PrimState m) -> m [Gen (PrimState m)]
splitGen n gen
  | n <= 0 = return []
  | otherwise = do
    seeds :: [V.Vector Word32] <- replicateM (n -1) $ uniformVector gen 256
    fmap (gen :) (mapM initialize seeds)

-- TODO: Splitmix. Remove or amend these functions as soon as split mix is used
-- and is available with the statistics package.

-- | Save a generator to a seed.
saveGen :: GenIO -> V.Vector Word32
saveGen = fromSeed . unsafePerformIO . save

-- | Load a generator from a seed.
loadGen :: V.Vector Word32 -> GenIO
loadGen = unsafePerformIO . restore . toSeed
