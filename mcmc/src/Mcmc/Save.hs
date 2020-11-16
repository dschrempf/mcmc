{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Save
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
-- Save and load an MCMC run. It is easy to save and restore the current state and
-- likelihood (or the trace), but it is not feasible to store all the proposals and so
-- on, so they have to be provided again when continuing a run.
module Mcmc.Save
  ( saveStatus,
    loadStatus,
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
-- TODO: Splitmix. Reproposal as soon as split mix is used and is available with the
-- statistics package.
import Mcmc.Item
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Chain hiding (save)
import Mcmc.Trace
import Mcmc.Verbosity
import Numeric.Log
import System.Directory
import System.IO.Unsafe (unsafePerformIO)
import System.Random.MWC
import Prelude hiding (cycle)

data Save a
  = Save
      -- Variables related to the chain.
      String -- Name.
      (Item a)
      Int -- Iteration.
      (Trace a)
      (Acceptance Int)
      (Maybe Int) -- Burn in.
      (Maybe Int) -- Auto tune.
      Int -- Iterations.
      Bool -- Force.
      (Maybe Int) -- Save.
      Verbosity
      (Vector Word32) -- Current seed.

      -- Variables related to the algorithm.
      [Maybe Double] -- Tuning parameters.

$(deriveJSON defaultOptions ''Save)

toSave :: Status a -> Save a
toSave (Status nm it i tr ac br at is f sv vb g _ _ _ _ _ c _) =
  Save
    nm
    it
    i
    tr'
    ac'
    br
    at
    is
    f
    sv
    vb
    g'
    ts
  where
    tr' = takeT (fromMaybe 0 sv) tr
    ac' = transformKeysA (ccProposals c) [0 ..] ac
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g' = fromSeed $ unsafePerformIO $ save g
    ts = [fmap tParam mt | mt <- map pTuner $ ccProposals c]

-- | Save a 'Status' to file.
--
-- Some important values have to be provided upon restoring the status. See
-- 'loadStatus'.
saveStatus :: ToJSON a => FilePath -> Status a -> IO ()
saveStatus fn s = BL.writeFile fn $ compress $ encode (toSave s)

-- fromSav prior lh cycle monitor save
fromSave ::
  (a -> Log Double) ->
  (a -> Log Double) ->
  Cycle a ->
  Monitor a ->
  Maybe (Cleaner a) ->
  Save a ->
  Status a
fromSave pr lh cc m cl (Save nm it i tr ac' br at is f sv vb g' ts) =
  Status
    nm
    it
    i
    tr
    ac
    br
    at
    is
    f
    sv
    vb
    g
    Nothing
    Nothing
    pr
    lh
    cl
    cc'
    m
  where
    ac = transformKeysA [0 ..] (ccProposals cc) ac'
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g = unsafePerformIO $ restore $ toSeed g'
    cc' = tuneCycle (M.mapMaybe id $ M.fromList $ zip (ccProposals cc) ts) cc

-- | Load a 'Status' from file.
--
-- Important information that cannot be saved and has to be provided again when
-- a chain is restored:
-- - prior function
-- - likelihood function
-- - cleaning function
-- - cycle
-- - monitor
--
-- To avoid incomplete continued runs, the @.mcmc@ file is removed after load.
loadStatus ::
  FromJSON a =>
  -- | Prior function.
  (a -> Log Double) ->
  -- | Likelihood function.
  (a -> Log Double) ->
  Cycle a ->
  Monitor a ->
  -- | Cleaner, if needed.
  Maybe (Cleaner a) ->
  -- | Path of status to load.
  FilePath ->
  IO (Status a)
loadStatus pr lh cc mn cl fn = do
  res <- eitherDecode . decompress <$> BL.readFile fn
  let s = case res of
        Left err -> error err
        Right sv -> fromSave pr lh cc mn cl sv
  -- Check if prior and likelihood matches.
  let Item x svp svl = item s
  -- Recompute and check the prior and likelihood for the last state because the
  -- functions may have changed. Of course, we cannot test for the same
  -- function, but having the same prior and likelihood at the last state is
  -- already a good indicator.
  when
    (pr x /= svp)
    (error "loadStatus: Provided prior function does not match the saved prior.")
  when
    (lh x /= svl)
    (error "loadStatus: Provided likelihood function does not match the saved likelihood.")
  removeFile fn
  return s
