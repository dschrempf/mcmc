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
    monitorBatchMeanInt,
    monitorBatchMeanIntF,
    monitorBatchMeanRealFloat,
    monitorBatchMeanRealFloatF,
    monitorBatchMeanRealFloatS,
    monitorBatchCustom,
  )
where

import Data.Text.Lazy.Builder (Builder)
import qualified Data.Text.Lazy.Builder.RealFloat as T
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
        mbpFunc :: [a] -> Builder
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

-- | Batch monitor integral parameters such as 'Int'. Print the mean with eight
-- decimal places (half precision).
monitorBatchMeanInt ::
  Integral a =>
  -- | Name of monitor.
  String ->
  MonitorParameterBatch a
monitorBatchMeanInt n = MonitorParameterBatch n (T.formatRealFloat T.Fixed (Just 8) . mean)
{-# SPECIALIZE monitorBatchMeanInt :: String -> MonitorParameterBatch Int #-}

-- | Batch monitor integral parameters such as 'Int'. Print the mean with full
-- precision.
monitorBatchMeanIntF ::
  Integral a =>
  -- | Name of monitor.
  String ->
  MonitorParameterBatch a
monitorBatchMeanIntF n = MonitorParameterBatch n (T.realFloat . mean)
{-# SPECIALIZE monitorBatchMeanIntF :: String -> MonitorParameterBatch Int #-}

-- | Batch monitor real float parameters such as 'Double' with eight decimal
-- places (half precision).
monitorBatchMeanRealFloat ::
  RealFloat a =>
  -- | Name of monitor.
  String ->
  MonitorParameterBatch a
monitorBatchMeanRealFloat n = MonitorParameterBatch n (T.formatRealFloat T.Fixed (Just 8) . mean)
{-# SPECIALIZE monitorBatchMeanRealFloat :: String -> MonitorParameterBatch Double #-}

-- | Batch monitor real float parameters such as 'Double' with full precision.
monitorBatchMeanRealFloatF ::
  RealFloat a =>
  -- | Name of monitor.
  String ->
  MonitorParameterBatch a
monitorBatchMeanRealFloatF n = MonitorParameterBatch n (T.realFloat . mean)
{-# SPECIALIZE monitorBatchMeanRealFloatF :: String -> MonitorParameterBatch Double #-}

-- | Batch monitor real float parameters such as 'Double' with scientific
-- notation and eight decimal places.
monitorBatchMeanRealFloatS ::
  RealFloat a =>
  -- | Name of monitor.
  String ->
  MonitorParameterBatch a
monitorBatchMeanRealFloatS n = MonitorParameterBatch n (T.formatRealFloat T.Exponent (Just 8) . mean)
{-# SPECIALIZE monitorBatchMeanRealFloatS :: String -> MonitorParameterBatch Double #-}

-- | Batch monitor parameters with custom lens and builder.
monitorBatchCustom ::
  -- | Name of monitor.
  String ->
  -- | Function to calculate the batch mean.
  ([a] -> a) ->
  -- | Custom builder.
  (a -> Builder) ->
  MonitorParameterBatch a
monitorBatchCustom n f b = MonitorParameterBatch n (b . f)
