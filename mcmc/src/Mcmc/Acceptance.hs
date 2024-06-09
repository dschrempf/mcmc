{-# LANGUAGE TemplateHaskell #-}

-- |
-- Module      :  Mcmc.Acceptance
-- Description :  Handle acceptance rates
-- Copyright   :  2021 Dominik Schrempf
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
    AcceptanceRates (..),
    Acceptance,
    Acceptances (fromAcceptances),
    emptyA,
    pushAccept,
    pushReject,
    ResetAcceptance (..),
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

-- | Proposals based on Hamiltonian dynamics use expected acceptance rates, not counts.
data AcceptanceRates = AcceptanceRates
  { totalAcceptanceRate :: !Double,
    nAcceptanceRates :: !Int
  }
  deriving (Show, Eq)

$(deriveJSON defaultOptions ''AcceptanceRates)

-- | Stored actual acceptance counts and maybe expected acceptance rates.
data Acceptance = A AcceptanceCounts (Maybe AcceptanceRates)
  deriving (Show, Eq)

$(deriveJSON defaultOptions ''Acceptance)

addAccept :: Maybe AcceptanceRates -> Acceptance -> Acceptance
addAccept mr' (A (AcceptanceCounts a r) mr) = A (AcceptanceCounts (a + 1) r) (addAcceptanceRates mr' mr)

addReject :: Maybe AcceptanceRates -> Acceptance -> Acceptance
addReject mr' (A (AcceptanceCounts a r) mr) = A (AcceptanceCounts a (r + 1)) (addAcceptanceRates mr' mr)

addAcceptanceRates :: Maybe AcceptanceRates -> Maybe AcceptanceRates -> Maybe AcceptanceRates
addAcceptanceRates Nothing Nothing = Nothing
addAcceptanceRates (Just r) Nothing = Just r
addAcceptanceRates Nothing (Just r) = Just r
addAcceptanceRates (Just (AcceptanceRates al rl)) (Just (AcceptanceRates ar rr)) =
  Just $ AcceptanceRates (al + ar) (rl + rr)

-- | For each key @k@, store the number of accepted and rejected proposals.
newtype Acceptances k = Acceptances {fromAcceptances :: M.Map k Acceptance}
  deriving (Eq, Show)

instance (ToJSONKey k) => ToJSON (Acceptances k) where
  toJSON (Acceptances m) = toJSON m
  toEncoding (Acceptances m) = toEncoding m

instance (Ord k, FromJSONKey k) => FromJSON (Acceptances k) where
  parseJSON v = Acceptances <$> parseJSON v

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
emptyA :: (Ord k) => [k] -> Acceptances k
emptyA ks = Acceptances $ M.fromList [(k, A noCounts Nothing) | k <- ks]
  where
    noCounts = AcceptanceCounts 0 0

-- | For key @k@, add an accept.
pushAccept :: (Ord k) => Maybe AcceptanceRates -> k -> Acceptances k -> Acceptances k
pushAccept mr k = Acceptances . M.adjust (addAccept mr) k . fromAcceptances

-- | For key @k@, add a reject.
pushReject :: (Ord k) => Maybe AcceptanceRates -> k -> Acceptances k -> Acceptances k
pushReject mr k = Acceptances . M.adjust (addReject mr) k . fromAcceptances

-- | Reset acceptance specification.
data ResetAcceptance
  = -- | Reset actual acceptance counts and expected acceptance rates.
    ResetEverything
  | -- | Only reset expected acceptance rates.
    ResetExpectedRatesOnly

-- | Reset acceptance counts.
resetA :: (Ord k) => ResetAcceptance -> Acceptances k -> Acceptances k
resetA ResetEverything = emptyA . M.keys . fromAcceptances
resetA ResetExpectedRatesOnly = Acceptances . M.map f . fromAcceptances
  where
    f (A cs _) = A cs Nothing

transformKeys :: (Ord k1, Ord k2) => [(k1, k2)] -> M.Map k1 v -> M.Map k2 v
transformKeys ks m = foldl' insrt M.empty ks
  where
    insrt m' (k1, k2) = M.insert k2 (m M.! k1) m'

-- | Transform keys using the given lists. Keys not provided will not be present
-- in the new 'Acceptance' variable.
transformKeysA :: (Ord k1, Ord k2) => [(k1, k2)] -> Acceptances k1 -> Acceptances k2
transformKeysA ks = Acceptances . transformKeys ks . fromAcceptances

-- | Compute acceptance counts, and actual and expected acceptances rates for a
-- specific proposal.
--
-- Return @Just (accepts, rejects, acceptance rate)@.
--
-- Return 'Nothing' if no proposals have been accepted or rejected (division by
-- zero).
acceptanceRate ::
  (Ord k) =>
  k ->
  Acceptances k ->
  -- | (nAccepts, nRejects, actualRate, expectedRate)
  (Int, Int, Maybe AcceptanceRate, Maybe AcceptanceRate)
acceptanceRate k a = case fromAcceptances a M.!? k of
  Just (A (AcceptanceCounts as rs) mrs) -> (as, rs, mar, mtr)
    where
      s = as + rs
      mar = if s <= 0 then Nothing else Just $ fromIntegral as / fromIntegral s
      mtr = case mrs of
        Nothing -> Nothing
        Just (AcceptanceRates xs n) -> Just $ xs / fromIntegral n
  Nothing -> error "acceptanceRate: Key not found in map."

-- | Compute actual acceptance rates for all proposals.
--
-- Set rate to 'Nothing' if no proposals have been accepted or rejected
-- (division by zero).
acceptanceRates :: Acceptances k -> M.Map k (Maybe AcceptanceRate)
acceptanceRates = M.map getRate . fromAcceptances
  where
    getRate (A (AcceptanceCounts as rs) _) =
      let s = as + rs
       in if s <= 0
            then Nothing
            else Just $ fromIntegral as / fromIntegral s
