-- |
-- Module      :  Mcmc.Proposal.SlideSpec
-- Description :  Unit tests for Mcmc.Proposal.SlideSpec
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May 19 12:07:43 2020.
module Mcmc.Proposal.SlideSpec
  ( spec,
  )
where

import Data.Maybe
import Mcmc.Proposal
import Mcmc.Proposal.Slide
import Test.Hspec
import Test.QuickCheck

prop_sym :: Eq b => (a -> a -> b) -> a -> a -> Bool
prop_sym f x y = f x y == f y x

slideSym :: Proposal Double
slideSym = slide "symmetric" 1 id 0 1.0 False

spec :: Spec
spec =
  describe "slide"
    $ it "has a symmetric proposal distribution if mean is 0"
    $ property
    $ prop_sym (fromJust $ mvDensity . mvSimple $ slideSym)
