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

module Mcmc.Trace
  ( Trace(..)
  , prependT
  )
where

import           Mcmc.Item

-- | A 'Trace' passes through a list of states with associated log-likelihoods
-- which are called 'Item's. New 'Item's are prepended, and the path of the
-- Markov chain is stored in reversed order.
newtype Trace a = Trace {fromTrace :: [Item a] }
  deriving (Show, Read)

instance Semigroup (Trace a) where
  (Trace l) <> (Trace r) = Trace (l <> r)

instance Monoid (Trace a) where
  mempty = Trace []

-- | Prepend an 'Item' to a 'Trace'.
prependT :: Item a -> Trace a -> Trace a
prependT x (Trace xs) = Trace (x : xs)
{-# INLINABLE prependT #-}

