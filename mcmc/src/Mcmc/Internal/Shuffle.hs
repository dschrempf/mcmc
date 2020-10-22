-- |
-- Module      :  Mcmc.Internal.Shuffle
-- Description :  Shuffle a list
-- Copyright   :  (c) Dominik Schrempf 2020
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
    shuffleN,
    grabble,
  )
where

import Control.Monad
import Control.Monad.ST
import Data.Vector (Vector)
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M
import System.Random.MWC
  ( GenIO,
    uniformR,
  )

-- | Shuffle a list.
shuffle :: [a] -> GenIO -> IO [a]
shuffle xs g = head <$> grabble xs 1 (length xs) g

-- | Shuffle a list @n@ times.
shuffleN :: [a] -> Int -> GenIO -> IO [[a]]
shuffleN xs n = grabble xs n (length xs)

-- -- Using System.Random.Shuffle. Speed is the same, so stay without additional dependency.
-- -- | Shuffle a list @n@ times.
-- shuffleN :: [a] -> Int -> GenIO -> IO [[a]]
-- shuffleN xs n g = replicateM n $ fmap (shuffle xs) (rseqM (length xs - 1) g)
--   where
--     rseqM :: Int -> GenIO -> IO [Int]
--     rseqM 0 _ = return []
--     rseqM i gen = liftM2 (:) (uniformR (0, i) gen) (rseqM (i - 1) gen)

-- | @grabble xs m n@ is /O(m*n')/, where @n' = min n (length xs)@. Choose @n'@
-- elements from @xs@, without replacement, and that @m@ times.
grabble :: [a] -> Int -> Int -> GenIO -> IO [[a]]
grabble xs m n gen = do
  swapss <- replicateM m $
    forM [0 .. min (l - 1) n] $ \i -> do
      j <- uniformR (i, l) gen
      return (i, j)
  return $ map (V.toList . V.take n . swapElems (V.fromList xs)) swapss
  where
    l = length xs - 1

swapElems :: Vector a -> [(Int, Int)] -> Vector a
swapElems xs swaps = runST $ do
  mxs <- V.unsafeThaw xs
  mapM_ (uncurry $ M.unsafeSwap mxs) swaps
  V.unsafeFreeze mxs