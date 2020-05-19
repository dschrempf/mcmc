{- |
Module      :  Statistics.Mcmc.Tools
Description :  Miscellaneous functions
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 12:31:40 2020.

-}

module Statistics.Mcmc.Tools
  ( acceptanceRatios
  , acceptanceRatio
  ) where

import qualified Data.Map.Strict as M
import Data.Map.Strict (Map)

ratio :: [Bool] -> Double
ratio xs = fromIntegral (length ts) / fromIntegral (length xs)
  where ts = filter (==True) xs

-- | Acceptance ratios averaged over the last @n@ iterations. If less than @n@
-- iterations are available, only those are used.
acceptanceRatios :: Int -> Map k [Bool] -> Map k Double
acceptanceRatios n = M.map (ratio . take n)

-- | Total acceptance ratio averaged over the last @n@ iterations. If less than
-- @n@ iterations are available, only those are used.
acceptanceRatio :: Int -> Map k [Bool] -> Double
acceptanceRatio n m = total / fromIntegral (M.size m)
  where total = sum (acceptanceRatios n m)
