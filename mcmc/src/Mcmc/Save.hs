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
-- Save and load chains. It is easy to save and restore the current state and
-- likelihood (or the trace), but it is not feasible to store all the proposals
-- and so on, so they have to be provided again when continuing a run.
module Mcmc.Save
  ( save,
    load,
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
import Mcmc.Chain
import Mcmc.Environment
import Mcmc.Item
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Trace
import Numeric.Log
import System.Directory
import System.IO.Unsafe (unsafePerformIO)
import qualified System.Random.MWC as MWC
import Prelude hiding (cycle)

data Save a
  = Save
      Environment
      String -- Name.
      (Item a)
      Int -- Iteration.
      (Trace a)
      (Acceptance Int)
      (Maybe Int) -- Burn in.
      (Maybe Int) -- Auto tune.
      Int -- Iterations.
      (Vector Word32) -- Current seed.
      [Maybe Double] -- Tuning parameters.

$(deriveJSON defaultOptions ''Save)

toSave :: Environment -> Chain a -> Save a
toSave env (Chain nm it i tr ac br at is g _ _ _ _ _ cc _) =
  Save
    env
    nm
    it
    i
    tr'
    ac'
    br
    at
    is
    g'
    ts
  where
    tr' =
      takeT
        ( case saveChain env of
            NoSave -> 0
            SaveN n -> n
        )
        tr
    ps = ccProposals cc
    ac' = transformKeysA ps [0 ..] ac
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g' = MWC.fromSeed $ unsafePerformIO $ MWC.save g
    ts = [fmap tParam mt | mt <- map pTuner ps]

-- | Save a chain with environment to file.
--
-- Some important values have to be provided upon restoring the status. See
-- 'load'.
save :: ToJSON a => Environment -> Chain a -> IO ()
save e c = BL.writeFile fn $ compress $ encode (toSave e c)
  where
    fn = name c ++ ".mcmc"

-- fromSav prior lh cycle monitor save
fromSave ::
  (a -> Log Double) ->
  (a -> Log Double) ->
  Cycle a ->
  Monitor a ->
  Maybe (Cleaner a) ->
  Save a ->
  (Environment, Chain a)
fromSave pr lh cc m cl (Save env nm it i tr ac' br at is g' ts) =
  ( env,
    Chain
      nm
      it
      i
      tr
      ac
      br
      at
      is
      g
      Nothing
      Nothing
      pr
      lh
      cl
      cc'
      m
  )
  where
    ac = transformKeysA [0 ..] (ccProposals cc) ac'
    -- TODO: Splitmix. Remove as soon as split mix is used and is available with
    -- the statistics package.
    g = unsafePerformIO $ MWC.restore $ MWC.toSeed g'
    cc' = tuneCycle (M.mapMaybe id $ M.fromList $ zip (ccProposals cc) ts) cc

-- | Load a chain with environment from file.
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
load ::
  FromJSON a =>
  -- | Prior function.
  (a -> Log Double) ->
  -- | Likelihood function.
  (a -> Log Double) ->
  Cycle a ->
  Monitor a ->
  -- | Cleaner, if needed.
  Maybe (Cleaner a) ->
  -- | Name of chain to load.
  FilePath ->
  IO (Environment, Chain a)
load pr lh cc mn cl nm = do
  res <- eitherDecode . decompress <$> BL.readFile fn
  let (e, c) = case res of
        Left err -> error err
        Right sv -> fromSave pr lh cc mn cl sv
  -- Check if prior and likelihood matches.
  let Item x svp svl = item c
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
  return (e, c)
  where fn = nm ++ ".mcmc"
