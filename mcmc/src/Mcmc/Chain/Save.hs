{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Chain.Save
-- Description :  Save the state of a Markov chain
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
  ( saveChainWith,
    loadChainWith,
  )
where

import Codec.Compression.GZip
import Control.Monad
import Data.Aeson
import Data.Aeson.TH
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List hiding (cycle)
import qualified Data.Map as M
import Data.Maybe
import Data.Vector.Unboxed (Vector)
import Data.Word
import Mcmc.Chain.Chain
import Mcmc.Chain.Item
import Mcmc.Chain.Trace
import Mcmc.Monitor
import Mcmc.Proposal
import System.Directory
import System.IO.Unsafe (unsafePerformIO)
import qualified System.Random.MWC as MWC
import Prelude hiding (cycle)

data SavedChain a = SavedChain
  { savedItem :: Item a,
    savedIteration :: Int,
    savedTrace :: Trace a,
    savedAcceptance :: Acceptance Int,
    savedSeed :: Vector Word32,
    savedTuningParameters :: [Maybe Double]
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SavedChain)

-- Save chain with trace of given maximum length.
toSavedChainWith :: Int -> Chain a -> SavedChain a
toSavedChainWith n (Chain it i tr ac g _ _ _ cc _) =
  SavedChain it i tr' ac' g' ts
  where
    tr' = takeT n tr
    ps = ccProposals cc
    ac' = transformKeysA ps [0 ..] ac
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g' = MWC.fromSeed $ unsafePerformIO $ MWC.save g
    ts = [fmap tParam mt | mt <- map pTuner ps]

-- | Save a Markov chain.
saveChainWith :: ToJSON a => Int -> FilePath -> Chain a -> IO ()
saveChainWith n fn c = BL.writeFile fn $ compress $ encode $ toSavedChainWith n c

-- Recompute and check the prior and likelihood for the last state because the
-- functions may have changed. Of course, we cannot test for the same
-- function, but having the same prior and likelihood at the last state is
-- already a good indicator.
loadSavedChain ::
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  SavedChain a ->
  Chain a
loadSavedChain pr lh cc mn (SavedChain it i tr ac' g' ts)
  | pr (state it) /= prior it =
    error "fromSave: Provided prior function does not match the saved prior."
  | lh (state it) /= likelihood it =
    error "fromSave: Provided likelihood function does not match the saved likelihood."
  | otherwise = Chain it i tr ac g i pr lh cc' mn
  where
    ac = transformKeysA [0 ..] (ccProposals cc) ac'
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g = unsafePerformIO $ MWC.restore $ MWC.toSeed g'
    cc' = tuneCycle (M.mapMaybe id $ M.fromList $ zip (ccProposals cc) ts) cc

-- | Load a Markov chain from a file.
loadChainWith ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  Monitor a ->
  FilePath ->
  IO (Chain a)
loadChainWith pr lh cc mn fn = do
  res <- eitherDecode . decompress <$> BL.readFile fn
  let c = case res of
        Left err -> error err
        Right c' -> loadSavedChain pr lh cc mn c'
  return c
