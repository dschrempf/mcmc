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
  )
where

import Control.Monad
import Control.Monad.Primitive
import qualified Data.Vector as V
import Data.Word
import System.Random.MWC

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
