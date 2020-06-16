{-# LANGUAGE OverloadedStrings #-}

{- |
Module      :  Mcmc.Save
Description :  Save the state of a Markov chain
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue Jun 16 10:18:54 2020.

-}

-- TODO: Continue an MCMC run. It is easy to save and restore the current state
-- and likelihood (or the trace), but it is not feasible to store all the moves
-- and so on, so they have to be provided again when continuing a run.

-- TODO: upon continuation: recompute and check the posterior for the last state
-- because the posterior function may have changed. Of course, we cannot test
-- for the same function, but having the same posterior at the last state is
-- already a good indicator.


module Mcmc.Save
  ()
where

import           Data.Aeson
import qualified Data.Vector.Unboxed           as V
import           Data.Vector.Unboxed            ( Vector )
import           Data.Word

import           Numeric.Log

import           Mcmc.Item
import           Mcmc.Monitor
import           Mcmc.Move
import           Mcmc.Status
import           Mcmc.Trace

import           System.Random.MWC

-- | Information about a Markov chain run, which can be (and is) stored on disk.
data Save a = Save (Item a) Int (Trace a) (Acceptance Int) (Vector Word32)

instance (ToJSON a) => ToJSON (Save a) where
  toJSON (Save s i t a g) =
    object ["s" .= s, "i" .= i, "t" .= t, "a" .= a, "g" .= g]
  toEncoding (Save s i t a g) =
    pairs ("s" .= s <> "i" .= i <> "t" .= t <> "a" .= a <> "g" .= g)

instance (FromJSON a) => FromJSON (Save a) where
  parseJSON = withObject "Save" $ \v ->
    Save <$> v .: "s" <*> v .: "i" <*> v .: "t" <*> v .: "a" <*> v .: "g"

toSave :: Status a -> Save a
toSave = undefined

-- fromSav logprior logllh cycle monitor save
fromSave :: (a -> Log Double) -> (a -> Log Double) -> Cycle a -> Monitor a -> Save a -> Status a
fromSave = undefined
