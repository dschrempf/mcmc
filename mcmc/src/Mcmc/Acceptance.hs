{-# LANGUAGE BangPatterns #-}

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
    Acceptance (fromAcceptance),
    emptyA,
    pushA,
    resetA,
    transformKeysA,
    acceptanceRate,
    acceptanceRates,
  )
where

import Control.DeepSeq
import Data.Aeson
import Data.Bifunctor
import Data.Foldable
import qualified Data.Map.Strict as M

-- | Acceptance rate.
type AcceptanceRate = Double

-- | For each key @k@, store the number of accepted and rejected proposals.
newtype Acceptance k = Acceptance {fromAcceptance :: M.Map k (Int, Int)}
  deriving (Eq, Read, Show)

instance ToJSONKey k => ToJSON (Acceptance k) where
  toJSON (Acceptance m) = toJSON m
  toEncoding (Acceptance m) = toEncoding m

instance (Ord k, FromJSONKey k) => FromJSON (Acceptance k) where
  parseJSON v = Acceptance <$> parseJSON v

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
emptyA :: Ord k => [k] -> Acceptance k
emptyA ks = Acceptance $ M.fromList [(k, (0, 0)) | k <- ks]

incrInt :: Int -> Int
incrInt x = x'
  where
    !x' = x + 1
{-# INLINE incrInt #-}

-- | For key @k@, prepend an accepted (True) or rejected (False) proposal.
pushA :: Ord k => k -> Bool -> Acceptance k -> Acceptance k
pushA k b =
  Acceptance
    . M.adjust
      (\(nA, nR) -> if b then (incrInt nA, nR) else (nA, incrInt nR))
      k
    . fromAcceptance
{-# INLINEABLE pushA #-}

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
  Just (0, 0) -> Nothing
  Just (as, rs) -> Just (as, rs, fromIntegral as / fromIntegral (as + rs))
  Nothing -> error "acceptanceRate: Key not found in map."

-- | Acceptance rates for all proposals.
--
-- Set rate to 'Nothing' if no proposals have been accepted or rejected
-- (division by zero).
acceptanceRates :: Acceptance k -> M.Map k (Maybe AcceptanceRate)
acceptanceRates =
  M.map
    ( \(as, rs) ->
        if as + rs == 0
          then Nothing
          else Just $ fromIntegral as / fromIntegral (as + rs)
    )
    . fromAcceptance
