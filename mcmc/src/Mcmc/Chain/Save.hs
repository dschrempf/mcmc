{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Jun 16 10:18:54 2020.
--
-- Save and load chains. It is easy to save and restore the current state and
-- likelihood (or the trace), but it is not feasible to store all the proposals
-- and so on, so they have to be provided again when continuing a run.

-- |
-- Module      :  Mcmc.Chain.Save
-- Description :  Save and load a Markov chain
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
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
import qualified Data.Stack.Circular as C
import qualified Data.Vector as VB
import qualified Data.Vector.Unboxed as VU
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
  { savedId :: Int,
    savedLink :: Link a,
    savedIteration :: Int,
    savedTrace :: C.Stack VB.Vector (Link a),
    savedAcceptance :: Acceptance Int,
    savedSeed :: VU.Vector Word32,
    savedTuningParameters :: [Maybe TuningParameter]
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SavedChain)

-- | Save a chain.
toSavedChain ::
  Chain a ->
  IO (SavedChain a)
toSavedChain (Chain ci it i tr ac g _ _ _ cc _) = do
  g' <- saveGen g
  tr' <- freezeT tr
  return $ SavedChain ci it i tr' ac' g' ts
  where
    ps = ccProposals cc
    ac' = transformKeysA ps [0 ..] ac
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
  IO (Chain a)
fromSavedChain pr lh cc mn (SavedChain ci it i tr ac' g' ts)
  | pr (state it) /= prior it =
    let msg =
          unlines
            [ "fromSave: Provided prior function does not match the saved prior.",
              "fromSave: Current prior:" <> show (prior it) <> ".",
              "fromSave: Given prior:" <> show (pr $ state it) <> "."
            ]
     in error msg
  | lh (state it) /= likelihood it =
    error "fromSave: Provided likelihood function does not match the saved likelihood."
  | otherwise = do
    g <- loadGen g'
    tr' <- thawT tr
    return $ Chain ci it i tr' ac g i pr lh cc' mn
  where
    ac = transformKeysA [0 ..] (ccProposals cc) ac'
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
