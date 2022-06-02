{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Acceptance
-- Description :  Handle acceptance rates
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Thu Jul  8 18:12:07 2021.
module Mcmc.Acceptance
  ( -- * Acceptance rates
    AcceptanceRate,
    AcceptanceCounts (..),
    Acceptance (fromAcceptance),
    emptyA,
    pushAccept,
    pushReject,
    pushAcceptanceCounts,
    resetA,
    transformKeysA,
    acceptanceRate,
    acceptanceRates,
  )
where

import Data.Aeson
import Data.Aeson.TH
import Data.Foldable
import qualified Data.Map.Strict as M

-- | Acceptance rate.
type AcceptanceRate = Double

-- | Number of accepted and rejected proposals.
data AcceptanceCounts = AcceptanceCounts
  { nAccepted :: !Int,
    nRejected :: !Int
  }
  deriving (Show, Eq, Ord)

$(deriveJSON defaultOptions ''AcceptanceCounts)

addAccept :: AcceptanceCounts -> AcceptanceCounts
addAccept (AcceptanceCounts a r) = AcceptanceCounts (a + 1) r

addReject :: AcceptanceCounts -> AcceptanceCounts
addReject (AcceptanceCounts a r) = AcceptanceCounts a (r + 1)

addAcceptanceCounts :: AcceptanceCounts -> AcceptanceCounts -> AcceptanceCounts
addAcceptanceCounts (AcceptanceCounts al rl) (AcceptanceCounts ar rr) =
  AcceptanceCounts (al + ar) (rl + rr)

-- | For each key @k@, store the number of accepted and rejected proposals.
newtype Acceptance k = Acceptance {fromAcceptance :: M.Map k AcceptanceCounts}
  deriving (Eq, Show)

instance ToJSONKey k => ToJSON (Acceptance k) where
  toJSON (Acceptance m) = toJSON m
  toEncoding (Acceptance m) = toEncoding m

instance (Ord k, FromJSONKey k) => FromJSON (Acceptance k) where
  parseJSON v = Acceptance <$> parseJSON v

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
emptyA :: Ord k => [k] -> Acceptance k
emptyA ks = Acceptance $ M.fromList [(k, AcceptanceCounts 0 0) | k <- ks]

-- | For key @k@, add an accept.
pushAccept :: Ord k => k -> Acceptance k -> Acceptance k
pushAccept k = Acceptance . M.adjust addAccept k . fromAcceptance

-- | For key @k@, add a reject.
pushReject :: Ord k => k -> Acceptance k -> Acceptance k
pushReject k = Acceptance . M.adjust addReject k . fromAcceptance

-- | For key @k@, add acceptance counts.
pushAcceptanceCounts :: Ord k => k -> AcceptanceCounts -> Acceptance k -> Acceptance k
pushAcceptanceCounts k c = Acceptance . M.adjust (addAcceptanceCounts c) k . fromAcceptance

-- | Reset acceptance storage.
resetA :: Ord k => Acceptance k -> Acceptance k
resetA = emptyA . M.keys . fromAcceptance

transformKeys :: (Ord k1, Ord k2) => [k1] -> [k2] -> M.Map k1 v -> M.Map k2 v
transformKeys ks1 ks2 m = foldl' insrt M.empty $ zip ks1 ks2
  where
    insrt m' (k1, k2) = M.insert k2 (m M.! k1) m'

-- | Transform keys using the given lists. Keys not provided will not be present
-- in the new 'Acceptance' variable.
transformKeysA :: (Ord k1, Ord k2) => [k1] -> [k2] -> Acceptance k1 -> Acceptance k2
transformKeysA ks1 ks2 = Acceptance . transformKeys ks1 ks2 . fromAcceptance

-- | Acceptance counts and rate for a specific proposal.
--
-- Return 'Nothing' if no proposals have been accepted or rejected (division by
-- zero).
acceptanceRate :: Ord k => k -> Acceptance k -> Maybe (Int, Int, AcceptanceRate)
acceptanceRate k a = case fromAcceptance a M.!? k of
  Just (AcceptanceCounts 0 0) -> Nothing
  Just (AcceptanceCounts as rs) -> Just (as, rs, fromIntegral as / fromIntegral (as + rs))
  Nothing -> error "acceptanceRate: Key not found in map."

-- | Acceptance rates for all proposals.
--
-- Set rate to 'Nothing' if no proposals have been accepted or rejected
-- (division by zero).
acceptanceRates :: Acceptance k -> M.Map k (Maybe AcceptanceRate)
acceptanceRates =
  M.map
    ( \(AcceptanceCounts as rs) ->
        if as + rs == 0
          then Nothing
          else Just $ fromIntegral as / fromIntegral (as + rs)
    )
    . fromAcceptance
