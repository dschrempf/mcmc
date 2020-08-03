{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Monitor.ParameterBatch
-- Description :  Monitor parameters
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 11:11:49 2020.
--
-- Batch mean monitors.
module Mcmc.Monitor.ParameterBatch
  ( -- * Batch parameter monitors
    MonitorParameterBatch (..),
    (@#),
    monitorBatchMean,
    monitorBatchMeanF,
    monitorBatchMeanE,
    monitorBatchCustom,
  )
where

import qualified Data.ByteString.Builder as BB
import qualified Data.Double.Conversion.ByteString as BC
import Lens.Micro

-- | Instruction about a parameter to monitor via batch means. Usually, the
-- monitored parameter is average over the batch size. However, arbitrary
-- functions performing more complicated analyses on the states in the batch can
-- be provided.
--
-- XXX: Batch monitors are slow at the moment because the monitored parameter
-- has to be extracted from the state for each iteration.
data MonitorParameterBatch a
  = MonitorParameterBatch
      { -- | Name of batch monitored parameter.
        mbpName :: String,
        -- | Instruction about how to extract the batch mean from the trace.
        mbpFunc :: [a] -> BB.Builder
      }

-- | Convert a batch parameter monitor from one data type to another using a
-- lens.
--
-- For example, to batch monitor a real float value being the first entry of a tuple:
--
-- @
-- mon = _1 @# monitorBatchMeanRealFloat
-- @
(@#) :: Lens' b a -> MonitorParameterBatch a -> MonitorParameterBatch b
(@#) l (MonitorParameterBatch n f) = MonitorParameterBatch n (f . map (^. l))

mean :: Real a => [a] -> Double
mean xs = realToFrac (sum xs) / fromIntegral (length xs)
{-# SPECIALIZE mean :: [Double] -> Double #-}
{-# SPECIALIZE mean :: [Int] -> Double #-}

-- | Batch monitor. Print the mean with eight decimal places (half precision).
monitorBatchMean ::
  Real a =>
  -- | Name.
  String ->
  MonitorParameterBatch a
monitorBatchMean n = MonitorParameterBatch n (BB.byteString . BC.toFixed 8 . mean)
{-# SPECIALIZE monitorBatchMean :: String -> MonitorParameterBatch Int #-}
{-# SPECIALIZE monitorBatchMean :: String -> MonitorParameterBatch Double #-}

-- | Batch monitor. Print the mean with full precision computing the shortest
-- string of digits that correctly represent the number.
monitorBatchMeanF ::
  Real a =>
  -- | Name.
  String ->
  MonitorParameterBatch a
monitorBatchMeanF n = MonitorParameterBatch n (BB.byteString . BC.toShortest . mean)
{-# SPECIALIZE monitorBatchMeanF :: String -> MonitorParameterBatch Int #-}
{-# SPECIALIZE monitorBatchMeanF :: String -> MonitorParameterBatch Double #-}

-- | Batch monitor real float parameters such as 'Double' with scientific
-- notation and eight decimal places.
monitorBatchMeanE ::
  Real a =>
  -- | Name.
  String ->
  MonitorParameterBatch a
monitorBatchMeanE n = MonitorParameterBatch n (BB.byteString . BC.toExponential 8 . mean)
{-# SPECIALIZE monitorBatchMeanE :: String -> MonitorParameterBatch Int #-}
{-# SPECIALIZE monitorBatchMeanE :: String -> MonitorParameterBatch Double #-}

-- | Batch monitor parameters with custom lens and builder.
monitorBatchCustom ::
  -- | Name.
  String ->
  -- | Function to calculate the batch mean.
  ([a] -> a) ->
  -- | Custom builder.
  (a -> BB.Builder) ->
  MonitorParameterBatch a
monitorBatchCustom n f b = MonitorParameterBatch n (b . f)
