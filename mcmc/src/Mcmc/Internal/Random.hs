{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      :  Mcmc.Internal.Random
-- Description :  Tools for random calculations
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Nov 25 07:14:52 2020.
module Mcmc.Internal.Random
  ( saveGen,
    loadGen,
  )
where

import Data.IORef
import Data.Word
import System.Random.Internal
import System.Random.SplitMix
import System.Random.Stateful

-- | Save a generator to a seed.
saveGen :: IOGenM StdGen -> IO (Word64, Word64)
saveGen (IOGenM r) = do
  (StdGen g) <- readIORef r
  pure $ unseedSMGen g

-- | Load a generator from a seed.
loadGen :: (Word64, Word64) -> IO (IOGenM StdGen)
loadGen s = newIOGenM $ StdGen $ seedSMGen' s
