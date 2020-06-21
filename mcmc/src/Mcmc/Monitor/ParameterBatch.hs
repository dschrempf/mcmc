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
  ( MonitorParameterBatch (..),
    monitorBatchMeanInt,
    monitorBatchMeanIntF,
    monitorBatchMeanRealFloat,
    monitorBatchMeanRealFloatF,
    monitorBatchMeanRealFloatS,
    monitorBatchCustom,
  )
where

import Data.Text.Lazy (Text)
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
        mbpName :: Text,
        -- | Instruction about how to extract the batch mean from the trace.
        mbpFunc :: [a] -> Builder
      }

mapL :: Lens' a b -> [a] -> [b]
mapL l = map (^. l)

mean :: Real a => [a] -> Double
mean xs = realToFrac (sum xs) / fromIntegral (length xs)
{-# SPECIALIZE mean :: [Double] -> Double #-}
{-# SPECIALIZE mean :: [Int] -> Double #-}

-- | Batch monitor integral parameters such as 'Int'. Print the mean with eight
-- decimal places (half precision).
monitorBatchMeanInt ::
  Integral b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameterBatch a
monitorBatchMeanInt n l =
  MonitorParameterBatch
    n
    (T.formatRealFloat T.Fixed (Just 8) . mean . mapL l)
{-# SPECIALIZE monitorBatchMeanInt :: Text -> Lens' a Int -> MonitorParameterBatch a #-}

-- | Batch monitor integral parameters such as 'Int'. Print the mean with full
-- precision.
monitorBatchMeanIntF ::
  Integral b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameterBatch a
monitorBatchMeanIntF n l =
  MonitorParameterBatch n (T.realFloat . mean . mapL l)
{-# SPECIALIZE monitorBatchMeanIntF :: Text -> Lens' a Int -> MonitorParameterBatch a #-}

-- | Batch monitor real float parameters such as 'Double' with eight decimal
-- places (half precision).
monitorBatchMeanRealFloat ::
  RealFloat b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameterBatch a
monitorBatchMeanRealFloat n l =
  MonitorParameterBatch
    n
    (T.formatRealFloat T.Fixed (Just 8) . mean . mapL l)
{-# SPECIALIZE monitorBatchMeanRealFloat :: Text -> Lens' a Double -> MonitorParameterBatch a #-}

-- | Batch monitor real float parameters such as 'Double' with full precision.
monitorBatchMeanRealFloatF ::
  RealFloat b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameterBatch a
monitorBatchMeanRealFloatF n l =
  MonitorParameterBatch n (T.realFloat . mean . mapL l)
{-# SPECIALIZE monitorBatchMeanRealFloatF :: Text -> Lens' a Double -> MonitorParameterBatch a #-}

-- | Batch monitor real float parameters such as 'Double' with scientific
-- notation and eight decimal places.
monitorBatchMeanRealFloatS ::
  RealFloat b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameterBatch a
monitorBatchMeanRealFloatS n l =
  MonitorParameterBatch
    n
    (T.formatRealFloat T.Exponent (Just 8) . mean . mapL l)
{-# SPECIALIZE monitorBatchMeanRealFloatS :: Text -> Lens' a Double -> MonitorParameterBatch a #-}

-- | Batch monitor parameters with custom lens and builder.
monitorBatchCustom ::
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  -- | Function to calculate the batch mean.
  ([b] -> b) ->
  -- | Custom builder.
  (b -> Builder) ->
  MonitorParameterBatch a
monitorBatchCustom n l f b =
  MonitorParameterBatch n (b . f . mapL l)
