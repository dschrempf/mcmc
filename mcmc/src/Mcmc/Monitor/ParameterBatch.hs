{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Monitor.ParameterBatch
-- Description :  Batch monitor parameters
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 11:11:49 2020.
--
-- A batch monitor prints summary statistics of a parameter collected over a
-- specific number of last iterations. The functions provided in this module
-- calculate the mean of the monitored parameter. However, custom batch monitors
-- can use more complex functions.
module Mcmc.Monitor.ParameterBatch
  ( -- * Batch parameter monitors
    MonitorParameterBatch (..),
    (>$<),
    monitorBatchMean,
    monitorBatchMeanF,
    monitorBatchMeanE,
  )
where

import qualified Data.ByteString.Builder as BB
import qualified Data.Double.Conversion.ByteString as BC
import Data.Functor.Contravariant
import qualified Data.Vector as VB

-- | Instruction about a parameter to monitor via batch means. Usually, the
-- monitored parameter is averaged over the batch size. However, arbitrary
-- functions performing more complicated analyses on the states in the batch can
-- be provided.
--
-- Convert a batch monitor from one data type to another with '(>$<)'.
--
-- For example, batch monitor the mean of the first entry of a tuple:
--
-- @
-- mon = fst >$< monitorBatchMean
-- @
--
-- Batch monitors may be slow because the monitored parameter has to be
-- extracted from the state for each iteration.
data MonitorParameterBatch a = MonitorParameterBatch
  { -- | Name of batch monitored parameter.
    mbpName :: String,
    -- | For a given batch, extract the summary statistics.
    mbpFunc :: VB.Vector a -> BB.Builder
  }

instance Contravariant MonitorParameterBatch where
  contramap f (MonitorParameterBatch n m) = MonitorParameterBatch n (m . VB.map f)

mean :: Real a => VB.Vector a -> Double
mean xs = realToFrac (VB.sum xs) / fromIntegral (VB.length xs)
{-# SPECIALIZE mean :: VB.Vector Double -> Double #-}
{-# SPECIALIZE mean :: VB.Vector Int -> Double #-}

-- | Batch mean monitor.
--
-- Print the mean with eight decimal places (half precision).
monitorBatchMean ::
  Real a =>
  -- | Name.
  String ->
  MonitorParameterBatch a
monitorBatchMean n = MonitorParameterBatch n (BB.byteString . BC.toFixed 8 . mean)
{-# SPECIALIZE monitorBatchMean :: String -> MonitorParameterBatch Int #-}
{-# SPECIALIZE monitorBatchMean :: String -> MonitorParameterBatch Double #-}

-- | Batch mean monitor.
--
-- Print the mean with full precision computing the shortest string of digits
-- that correctly represent the number.
monitorBatchMeanF ::
  Real a =>
  -- | Name.
  String ->
  MonitorParameterBatch a
monitorBatchMeanF n = MonitorParameterBatch n (BB.byteString . BC.toShortest . mean)
{-# SPECIALIZE monitorBatchMeanF :: String -> MonitorParameterBatch Int #-}
{-# SPECIALIZE monitorBatchMeanF :: String -> MonitorParameterBatch Double #-}

-- | Batch mean monitor.
--
-- Print the real float parameters such as 'Double' with scientific notation and
-- eight decimal places.
monitorBatchMeanE ::
  Real a =>
  -- | Name.
  String ->
  MonitorParameterBatch a
monitorBatchMeanE n = MonitorParameterBatch n (BB.byteString . BC.toExponential 8 . mean)
{-# SPECIALIZE monitorBatchMeanE :: String -> MonitorParameterBatch Int #-}
{-# SPECIALIZE monitorBatchMeanE :: String -> MonitorParameterBatch Double #-}
