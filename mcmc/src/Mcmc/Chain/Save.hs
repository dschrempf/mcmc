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
  ( save,
    loadWith,
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
-- TODO: Splitmix. Remove as soon as split mix is used and is available with the
-- statistics package.
import Data.Word
import Mcmc.Chain.Chain
import Mcmc.Chain.Item
import Mcmc.Chain.Trace
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import System.Directory
import System.IO.Unsafe (unsafePerformIO)
import qualified System.Random.MWC as MWC
import Prelude hiding (cycle)

data SavedChain a = SavedChain
  { savedChainName :: String,
    savedItem :: Item a,
    savedIteration :: Int,
    savedTrace :: Trace a,
    savedAcceptance :: Acceptance Int,
    savedSeed :: Vector Word32,
    savedTuningParameters :: [Maybe Double]
  }
  deriving (Eq, Read, Show)

$(deriveJSON defaultOptions ''SavedChain)

data Saved a = Saved
  { savedSettings :: Settings,
    savedChains :: [SavedChain a]
  }
  deriving (Eq, Show)

$(deriveJSON defaultOptions ''Saved)

-- Save chain with trace of given maximum length.
saveChainWith :: Int -> Chain a -> SavedChain a
saveChainWith n (Chain nm it i tr ac g _ _ _ cc _) =
  SavedChain nm it i tr' ac' g' ts
  where
    tr' = takeT n tr
    ps = ccProposals cc
    ac' = transformKeysA ps [0 ..] ac
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g' = MWC.fromSeed $ unsafePerformIO $ MWC.save g
    ts = [fmap tParam mt | mt <- map pTuner ps]

-- | Save an MCMC run to file.
--
-- Some important values have to be provided upon 'load'.
save :: ToJSON a => Settings -> [Chain a] -> IO ()
save e cs = BL.writeFile fn $ compress $ encode $ Saved e $ map (saveChainWith nTr) cs
  where
    fn = name e ++ ".mcmc"
    nTr = case saveMode e of
      NoSave -> 0
      SaveWithTrace n -> n

-- loadChain prior lh cycle monitor save
--
-- Recompute and check the prior and likelihood for the last state because the
-- functions may have changed. Of course, we cannot test for the same
-- function, but having the same prior and likelihood at the last state is
-- already a good indicator.
loadChain ::
  PriorFunction a ->
  LikelihoodFunction a ->
  -- CleaningFunction a ->
  Cycle a ->
  Monitor a ->
  SavedChain a ->
  Chain a
loadChain pr lh cc mn (SavedChain nm it i tr ac' g' ts)
  | fromPriorFunction pr (state it) /= prior it =
    error "fromSave: Provided prior function does not match the saved prior."
  | fromLikelihoodFunction lh (state it) /= likelihood it =
    error "fromSave: Provided likelihood function does not match the saved likelihood."
  | otherwise = Chain nm it i tr ac g i pr lh cc' mn
  where
    ac = transformKeysA [0 ..] (ccProposals cc) ac'
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g = unsafePerformIO $ MWC.restore $ MWC.toSeed g'
    cc' = tuneCycle (M.mapMaybe id $ M.fromList $ zip (ccProposals cc) ts) cc

-- TODO: REFACTOR. Different prior functions need to be given for different
-- chains in the list. Better define save and load with the 'Algorithm'.

-- | Load an MCMC run from file.
--
-- To avoid incomplete continued runs, the @.mcmc@ file is removed after load.
loadWith ::
  FromJSON a =>
  PriorFunction a ->
  LikelihoodFunction a ->
  -- CleaningFunction a ->
  Cycle a ->
  Monitor a ->
  -- | Name of chain to load.
  FilePath ->
  IO (Settings, [Chain a])
loadWith pr lh cc mn nm = do
  res <- eitherDecode . decompress <$> BL.readFile fn
  let (e, cs) = case res of
        Left err -> error err
        Right (Saved env cs') -> (env, map (loadChain pr lh cc mn) cs')
  removeFile fn
  return (e, cs)
  where
    fn = nm ++ ".mcmc"
