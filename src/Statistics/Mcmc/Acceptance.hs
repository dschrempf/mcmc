{- |
Module      :  Statistics.Mcmc.Acceptance
Description :  Miscellaneous functions
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May  6 12:31:40 2020.

-}

-- XXX: I guess using mutable vectors or something similar would be much faster.

module Statistics.Mcmc.Acceptance
  ( Acceptance
  , empty
  , prependA
  , acceptanceRatios
  , acceptanceRatio
  ) where

import qualified Data.Map.Strict as M
import Data.Map.Strict (Map)

-- | For each key @k@, store the list of accepted (True) and rejected (False)
-- proposals. For reasons of efficiency, the lists are stored in reverse order;
-- latest first.
newtype Acceptance k = Acceptance { fromAcceptance :: Map k [Bool] }

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
empty :: Ord k => [k] -> Acceptance k
empty ks = Acceptance $ M.fromList [ (k, []) | k <- ks ]

-- | For key @k@, prepend an accepted (True) or rejected (False) proposal.
prependA :: (Ord k, Show k) => k -> Bool -> Acceptance k -> Acceptance k
-- XXX: Unsafe; faster.
prependA k v (Acceptance m) = Acceptance $ M.adjust (v:) k m
-- -- XXX: Safe; slower.
-- prependA k v (Acceptance m) | k `M.member` m = Acceptance $ M.adjust (v:) k m
--                             | otherwise = error msg
--   where msg = "prependA: Can not add acceptance value for key: " <> show k <> "."
{-# INLINEABLE #-}

ratio :: [Bool] -> Double
ratio xs = fromIntegral (length ts) / fromIntegral (length xs)
  where ts = filter (==True) xs

-- | Acceptance ratios, averaged over the last @n@ iterations. If less than @n@
-- iterations are available, only those are used.
acceptanceRatios :: Int -> Acceptance k -> Map k Double
acceptanceRatios n (Acceptance m)= M.map (ratio . take n) m

-- | Total acceptance ratio averaged over the last @n@ iterations. If less than
-- @n@ iterations are available, only those are used.
acceptanceRatio :: Int -> Acceptance k -> Double
acceptanceRatio n a = total / fromIntegral (M.size $ fromAcceptance a)
  where total = sum (acceptanceRatios n a)

