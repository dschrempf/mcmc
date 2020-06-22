{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- |
-- Module      :  Mcmc.Item
-- Description :  Links of Markov chains
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 09:10:27 2020.
module Mcmc.Item
  ( Item (..),
  )
where

import Data.Aeson
import Numeric.Log

instance ToJSON a => ToJSON (Log a) where
  toJSON (Exp x) = toJSON x
  toEncoding (Exp x) = toEncoding x

instance FromJSON a => FromJSON (Log a) where
  parseJSON v = Exp <$> parseJSON v

-- | An 'Item', or link of the Markov chain. For reasons of computational
-- efficiency, each state is associated with the corresponding prior and
-- likelihood.
data Item a
  = Item
      { -- | The current state in the state space @a@.
        state :: a,
        -- | The current prior.
        prior :: Log Double,
        -- | The current likelihood.
        likelihood :: Log Double
      }
  deriving (Eq, Ord, Show, Read)

instance ToJSON a => ToJSON (Item a) where
  toJSON (Item x p l) = object ["s" .= x, "p" .= p, "l" .= l]
  toEncoding (Item x p l) = pairs ("s" .= x <> "p" .= p <> "l" .= l)

instance FromJSON a => FromJSON (Item a) where
  parseJSON = withObject "Item" $
    \v ->
      Item
        <$> v .: "s"
        <*> v .: "p"
        <*> v .: "l"
