{- |
Module      :  Statistics.Mcmc.Tools.Shuffle
Description :  Shuffle a list
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May 20 14:37:09 2020.

From https://wiki.haskell.org/Random_shuffle.

-}

module Statistics.Mcmc.Tools.Shuffle
  ( shuffle
  , shuffleN
  , grabble
  ) where

import           Control.Monad
import           Data.Array                     ( elems )
import           Data.Array.ST                  ( newListArray
                                                , readArray
                                                , runSTArray
                                                , writeArray
                                                )
import           System.Random.MWC

-- | Shuffle a list.
shuffle :: [a] -> GenIO -> IO [a]
shuffle xs g = head <$> grabble xs 1 (length xs) g

-- | Shuffle a list @n@ times.
shuffleN :: [a] -> Int -> GenIO -> IO [[a]]
shuffleN xs n = grabble xs n (length xs)

-- | @grabble xs m n@ is /O(m*n')/, where @n' = min n (length xs)@. Choose @n'@
-- elements from @xs@, without replacement, and that @m@ times.
grabble :: [a] -> Int -> Int -> GenIO -> IO [[a]]
grabble xs m n gen = do
  swapss <- replicateM m $ forM [0 .. min (maxIx - 1) n] $ \i -> do
    j <- uniformR (i, maxIx) gen
    return (i, j)
  return $ map (take n . swapElems xs) swapss
  where maxIx = length xs - 1

swapElems :: [a] -> [(Int, Int)] -> [a]
swapElems xs swaps = elems $ runSTArray
  (do
    arr <- newListArray (0, maxIx) xs
    mapM_ (swap arr) swaps
    return arr
  )
 where
  maxIx = length xs - 1
  swap arr (i, j) = do
    vi <- readArray arr i
    vj <- readArray arr j
    writeArray arr i vj
    writeArray arr j vi
