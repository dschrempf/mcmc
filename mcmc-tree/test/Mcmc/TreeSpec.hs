{-# LANGUAGE OverloadedStrings #-}

-- |
--   Module      :  Mcmc.TreeSpec
--   Description :  Unit tests for Mcmc.Tree
--   Copyright   :  (c) Dominik Schrempf, 2020
--   License     :  GPL-3.0-or-later
--
--   Maintainer  :  dominik.schrempf@gmail.com
--   Stability   :  unstable
--   Portability :  portable
--
-- Creation date: Fri Jul  3 10:19:55 2020.
module Mcmc.TreeSpec
  ( spec,
  )
where

import Test.Hspec
import Mcmc.Tree

fnE :: FilePath
fnE = "data/easy.tree"

fn :: FilePath
fn = "data/plants.treelist"

spec :: Spec
spec =
  describe "manyNewick" $ do
    it "reads an easy tree" $ do
      trs <- manyNewick fnE
      -- print $ head trs
      length trs `shouldBe` 1
    it "reads many trees" $ do
      trs <- manyNewick fn
      -- print $ head trs
      length trs `shouldBe` 2000
