{- |
Module      :  Statistics.Mcmc.MovesSpec
Description :  Unit tests for Statistics.Mcmc.MovesSpec
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May 19 12:07:43 2020.

-}

module Statistics.Mcmc.MovesSpec
  ( spec
  ) where

import           Test.Hspec
import           Test.QuickCheck


import Statistics.Mcmc.Moves
import Statistics.Mcmc.Types

prop_sym :: Eq b => (a -> a -> b) -> a -> a -> Bool
prop_sym f x y = f x y == f y x

slideSym :: Move Double
slideSym = slideDouble "symmetric" 0 1.0

spec :: Spec
spec = do
  describe "slide" $ do
    it "has a symmetric proposal distribution if mean is 0" $ property $
      prop_sym (mvLogDensity slideSym)
