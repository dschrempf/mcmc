{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Chain.Save
-- Description :  Save and load a Markov chain
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Jun 16 10:18:54 2020.
--
-- Save and load chains. It is easy to save and restore the current state and
-- likelihood (or the trace), but it is not feasible to store all the proposals
-- and so on, so they have to be provided again when continuing a run.
module Mcmc.Chain.Save
  ( SavedChain (..),
    toSavedChain,
    fromSavedChain,
  )
where

import Control.Monad
import Data.Aeson
import Data.Aeson.TH
import Data.List hiding (cycle)
import qualified Data.Map as M
import Data.Maybe
import Data.Vector.Unboxed (Vector)
import Data.Word
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Trace
import Mcmc.Internal.Random
import Mcmc.Monitor
import Mcmc.Proposal
import Prelude hiding (cycle)

-- | Storable values of a Markov chain.
--
-- See 'toSavedChain'.
data SavedChain a = SavedChain
  { savedLink :: Link a,
    savedIteration :: Int,
    savedTrace :: Trace a,
    savedAcceptance :: Acceptance Int,
    savedSeed :: Vector Word32,
    savedTuningParameters :: [Maybe Double]
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SavedChain)

-- | Save a chain.
toSavedChain ::
  -- | Maximum length of trace.
  Int ->
  Chain a ->
  SavedChain a
toSavedChain n (Chain it i tr ac g _ _ _ cc _) =
  SavedChain it i tr' ac' g' ts
  where
    tr' = takeT n tr
    ps = ccProposals cc
    ac' = transformKeysA ps [0 ..] ac
    g' = saveGen g
    ts = [fmap tParam mt | mt <- map pTuner ps]

-- | Load a saved chain.
--
-- Recompute and check the prior and likelihood for the last state because the
-- functions may have changed. Of course, we cannot test for the same function,
-- but having the same prior and likelihood at the last state is already a good
-- indicator.
fromSavedChain ::
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  SavedChain a ->
  Chain a
fromSavedChain pr lh cc mn (SavedChain it i tr ac' g' ts)
  | pr (state it) /= prior it =
    error "fromSave: Provided prior function does not match the saved prior."
  | lh (state it) /= likelihood it =
    error "fromSave: Provided likelihood function does not match the saved likelihood."
  | otherwise = Chain it i tr ac g i pr lh cc' mn
  where
    ac = transformKeysA [0 ..] (ccProposals cc) ac'
    g = loadGen g'
    getTuningF mt = case mt of
      Nothing -> const 1.0
      Just t -> const t
    cc' =
      tuneCycle
        ( M.map getTuningF $
            M.fromList $
              zip (ccProposals cc) ts
        )
        cc
