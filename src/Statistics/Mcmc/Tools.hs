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

import Statistics.Mcmc.Types

ratio :: [Bool] -> Double
ratio xs = fromIntegral (length ts) / fromIntegral (length xs)
  where ts = filter (==True) xs

-- TODO: Instead of passing a status, pass an (hypothetical) Acceptance type.
acceptanceRatios :: Status a -> [Double]
acceptanceRatios s = map ratio as
  where as = acceptance s

acceptanceRatio :: Status a -> Double
acceptanceRatio s = sum rs / fromIntegral (length rs)
  where rs = acceptanceRatios s
