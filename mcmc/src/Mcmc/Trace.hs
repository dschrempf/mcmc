-- |
-- Module      :  Mcmc.Trace
-- Description :  Trace of a Markov chain
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 09:11:25 2020.
module Mcmc.Trace
  ( Trace,
    singletonT,
    pushT,
    headT,
    takeT,
  )
where

import Data.Aeson
import Mcmc.Item

-- | A 'Trace' passes through a list of states with associated likelihoods which
-- are called 'Item's. New 'Item's are prepended, and the path of the Markov
-- chain is stored in reversed order.
newtype Trace a = Trace {fromTrace :: [Item a]}
  deriving (Show, Read, Eq)

instance Semigroup (Trace a) where
  (Trace l) <> (Trace r) = Trace (l <> r)

instance Monoid (Trace a) where
  mempty = Trace []

instance ToJSON a => ToJSON (Trace a) where
  toJSON (Trace xs) = toJSON xs
  toEncoding (Trace xs) = toEncoding xs

instance FromJSON a => FromJSON (Trace a) where
  parseJSON v = Trace <$> parseJSONList v

-- | The empty trace.
singletonT :: Item a -> Trace a
singletonT i = Trace [i]

-- | Prepend an 'Item' to a 'Trace'.
pushT :: Item a -> Trace a -> Trace a
pushT x = Trace . (:) x . fromTrace
{-# INLINEABLE pushT #-}

-- | Get the most recent item of the trace.
headT :: Trace a -> Item a
headT = head . fromTrace

-- | Get the N most recent items of the trace.
takeT :: Int -> Trace a -> [Item a]
takeT n = take n . fromTrace
