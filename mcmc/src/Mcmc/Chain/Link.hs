{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Chain.Link
-- Description :  The state combined with auxiliary variables
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 09:10:27 2020.
module Mcmc.Chain.Link
  ( Link (..),
  )
where

import Data.Aeson
import Data.Aeson.Types
import Mcmc.Likelihood
import Mcmc.Prior
import Numeric.Log

-- | Link of a Markov chain. For reasons of computational efficiency, each state
-- is associated with the corresponding prior and likelihood.
data Link a = Link
  { -- | The current state in the state space @a@.
    state :: a,
    -- | The current prior value.
    prior :: Prior,
    -- | The current likelihood value.
    likelihood :: Likelihood
  }
  deriving (Eq, Ord, Show, Read)

instance (ToJSON a) => ToJSON (Link a) where
  toJSON (Link x (Exp p) (Exp l)) = object ["s" .= x, "p" .= p, "l" .= l]
  toEncoding (Link x (Exp p) (Exp l)) = pairs ("s" .= x <> "p" .= p <> "l" .= l)

link :: (FromJSON a) => Object -> Parser (Link a)
link v = do
  s <- v .: "s"
  p <- v .: "p"
  l <- v .: "l"
  return $ Link s (Exp p) (Exp l)

instance (FromJSON a) => FromJSON (Link a) where
  parseJSON = withObject "Link" link
