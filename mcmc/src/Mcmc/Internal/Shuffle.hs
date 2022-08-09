-- |
-- Module      :  Mcmc.Internal.Shuffle
-- Description :  Shuffle a list
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 14:37:09 2020.
--
-- From https://wiki.haskell.org/Random_shuffle.
module Mcmc.Internal.Shuffle
  ( shuffle,
  )
where

import Control.Monad
import Control.Monad.ST
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M
import System.Random.Stateful

-- TODO: Remove shuffle, use System.Random.MWC.Distributions.uniformShuffle and
-- vectors.

-- | Shuffle a vector.
shuffle :: StatefulGen g m => [a] -> g -> m [a]
shuffle xs = grabble xs (length xs)

-- @grabble xs m n@ is /O(m*n')/, where @n' = min n (length xs)@. Choose @n'@
-- elements from @xs@, without replacement, and that @m@ times.
grabble :: StatefulGen g m => [a] -> Int -> g -> m [a]
grabble xs m gen = do
  swaps <- forM [0 .. min (l - 1) m] $ \i -> do
    j <- uniformRM (i, l) gen
    return (i, j)
  return $ (V.toList . V.take m . swapElems (V.fromList xs)) swaps
  where
    l = length xs - 1

swapElems :: V.Vector a -> [(Int, Int)] -> V.Vector a
swapElems xs swaps = runST $ do
  mxs <- V.unsafeThaw xs
  mapM_ (uncurry $ M.unsafeSwap mxs) swaps
  V.unsafeFreeze mxs
