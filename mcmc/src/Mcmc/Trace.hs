{- |
Module      :  Mcmc.Trace
Description :  Trace of a Markov chain
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May 20 09:11:25 2020.

-}

-- TODO: Possibly limit the length of the trace to the maximum batch size.
--
-- This could be achieved with a special data structure storing a fixed sized
-- vector together with the current index. However, I couldn't find such a
-- structure on Hackage.

module Mcmc.Trace
  ( Trace (fromTrace)
  , singleton
  , prependT
  , pop
  , takeN
  )
where

import           Data.Aeson

import           Mcmc.Item

-- | A 'Trace' passes through a list of states with associated log-likelihoods
-- which are called 'Item's. New 'Item's are prepended, and the path of the
-- Markov chain is stored in reversed order.
newtype Trace a = Trace {fromTrace :: [Item a] }
  deriving (Show, Read, Eq)

instance Semigroup (Trace a) where
  (Trace l) <> (Trace r) = Trace (l <> r)

instance Monoid (Trace a) where
  mempty = Trace []

instance ToJSON a => ToJSON (Trace a) where
  toJSON (Trace xs) = toJSON xs
  toEncoding (Trace xs) = toEncoding xs

instance FromJSON a => FromJSON (Trace a)  where
  parseJSON v = Trace <$> parseJSONList v

-- $(deriveJSON defaultOptions 'Trace)

-- | The empty trace.
singleton :: Item a -> Trace a
singleton i = Trace [i]

-- | Prepend an 'Item' to a 'Trace'.
prependT :: Item a -> Trace a -> Trace a
prependT x (Trace xs) = Trace (x : xs)
{-# INLINABLE prependT #-}

-- | Get the most recent item of the trace.
pop :: Trace a -> Item a
pop = head . fromTrace

-- | Get the N most recent item of the trace.
takeN :: Int -> Trace a -> Trace a
takeN n = Trace . take n . fromTrace
