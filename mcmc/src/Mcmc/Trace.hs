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

-- XXX: Trace provided as an abstract data type, because the implementation
-- should change from a list to a fixed sized vector with access index, see
-- above.

module Mcmc.Trace
  ( Trace
  , singletonT
  , prependT
  , headT
  , takeT
  , states
  , logPriors
  , logLikelihoods
  )
where

import           Data.Aeson
import           Numeric.Log

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
singletonT :: Item a -> Trace a
singletonT i = Trace [i]

-- | Prepend an 'Item' to a 'Trace'.
prependT :: Item a -> Trace a -> Trace a
prependT x (Trace xs) = Trace (x : xs)
{-# INLINABLE prependT #-}

-- | Get the most recent item of the trace.
headT :: Trace a -> Item a
headT = head . fromTrace

-- | Get the N most recent items of the trace.
takeT :: Int -> Trace a -> Trace a
takeT n = Trace . take n . fromTrace

-- | Get the actual states.
states :: Trace a -> [a]
states = map state . fromTrace

-- | Get the log priors.
logPriors :: Trace a -> [Log Double]
logPriors = map logPrior . fromTrace

-- | Get the log likelihoods.
logLikelihoods :: Trace a -> [Log Double]
logLikelihoods = map logPrior . fromTrace
