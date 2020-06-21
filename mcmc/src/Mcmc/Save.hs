{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

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
-- likelihood (or the trace), but it is not feasible to store all the moves and so
-- on, so they have to be provided again when continuing a run.
module Mcmc.Save
  ( saveStatus
  , loadStatus
  )
where

import           Control.Monad
import           Data.Aeson
import           Data.Aeson.TH
import           Data.List               hiding ( cycle )
import           Data.Map                       ( Map )
import qualified Data.Map                      as M
import           Data.Time.Clock
import           Data.Vector.Unboxed            ( Vector )
import           Data.Word
-- TODO: Remove as soon as split mix is used and is available with the
-- statistics package.
import           Mcmc.Item
import           Mcmc.Monitor
import           Mcmc.Move
import           Mcmc.Status
import           Mcmc.Trace
import           Numeric.Log
import           Prelude                 hiding ( cycle )
import           System.IO.Unsafe               ( unsafePerformIO )
import           System.Random.MWC

-- | Information about a Markov chain run, which can be (and is) stored on disk.
-- data Save a =
--   Save { sName :: String
--        , sItem :: Item a
--        , sIteration :: Int
--        , sTrace :: Trace a
--        , sAcceptance :: Acceptance Int
--        , sBurnInIterations :: Maybe Int
--        , sAutoTuningPeriod :: Maybe Int
--        , sIterations :: Int
--        , sStarttime :: Maybe UTCTime
--        , sSavetime :: Maybe UTCTime
--        , sGenerator :: Vector Word32
--        }
data Save a =
  Save String
       (Item a)
       Int
       (Trace a)
       (Acceptance Int)
       (Maybe Int)
       (Maybe Int)
       Int
       (Maybe UTCTime)
       (Maybe UTCTime)
       (Vector Word32)

$(deriveJSON defaultOptions 'Save)

mapKeys :: (Ord k1, Ord k2) => [(k1, k2)] -> Map k1 v -> Map k2 v
mapKeys xs m = foldl' addValue M.empty xs
  where addValue m' (k1, k2) = M.insert k2 (m M.! k1) m'

toSave :: Status a -> Save a
toSave (Status nm it i tr ac br at is st sv g _ _ c _) = Save nm
                                                              it
                                                              i
                                                              tr
                                                              ac'
                                                              br
                                                              at
                                                              is
                                                              st
                                                              sv
                                                              g'
 where
  mvs = flip zip [0 ..] $ fromCycle c

  ac' = Acceptance $ mapKeys mvs (fromAcceptance ac)

  -- TODO: Remove as soon as split mix is used and is available with the
  -- statistics package.
  g'  = fromSeed $ unsafePerformIO $ save g

-- | Save a 'Status' to file.
--
-- Saved information:
-- - state
-- - iteration
-- - trace
-- - acceptance ratios
-- - generator
--
-- Important information that cannot be saved and has to be provided again when
-- a chain is restored:
-- - log prior function
-- - log likelihood function
-- - cycle
-- - monitor
saveStatus :: ToJSON a => FilePath -> Status a -> IO ()
saveStatus fn s = encodeFile fn (toSave s)

-- fromSav logprior logllh cycle monitor save
fromSave
  :: (a -> Log Double)
  -> (a -> Log Double)
  -> Cycle a
  -> Monitor a
  -> Save a
  -> Status a
fromSave p l c m (Save nm it i tr ac br at is st sv g) = Status nm
                                                                it
                                                                i
                                                                tr
                                                                ac'
                                                                br
                                                                at
                                                                is
                                                                st
                                                                sv
                                                                g'
                                                                p
                                                                l
                                                                c
                                                                m
 where
  mvs = zip [0 ..] $ fromCycle c

  ac' = Acceptance $ mapKeys mvs (fromAcceptance ac)

  -- TODO: Remove as soon as split mix is used and is available with the
  -- statistics package.
  g'  = unsafePerformIO $ restore $ toSeed g

-- | Load a 'Status' from file.
-- Important information that cannot be saved and has to be provided again when
-- a chain is restored:
-- - log prior function
-- - log likelihood function
-- - cycle
-- - monitor
loadStatus
  :: FromJSON a
  => (a -> Log Double)
  -> (a -> Log Double)
  -> Cycle a
  -> Monitor a
  -> FilePath
  -> IO (Status a)
loadStatus p l c m fn = do
  res <- eitherDecodeFileStrict' fn
  let s = case res of
        Left  err -> error err
        Right sv  -> fromSave p l c m sv
  -- Check if log prior and log likelihood matches.
  let Item x svp svl = item s
  -- Recompute and check the prior and likelihood for the last state because the
  -- functions may have changed. Of course, we cannot test for the same
  -- function, but having the same prior and likelihood at the last state is
  -- already a good indicator.
  when
    (p x /= svp)
    (error
      "loadStatus: Provided log prior function does not match the saved log prior."
    )
  when
    (l x /= svl)
    (error
      "loadStatus: Provided log likelihood function does not match the saved log likelihood."
    )
  return s