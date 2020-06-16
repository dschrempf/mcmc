{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Mcmc.Monitor.Batch
Description :  Monitor parameters
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May 29 11:11:49 2020.

Batch mean monitors.

-}

module Mcmc.Monitor.Batch
  ( MonitorParameter(..)
  , monitorBatchInt
  , monitorBatchRealFloat
  , monitorBatchRealFloatF
  , monitorBatchRealFloatS
  )
where

-- TODO!

import qualified Data.Text.Lazy.Builder.Int    as T
import qualified Data.Text.Lazy.Builder.RealFloat
                                               as T
import           Data.Text.Lazy                 ( Text )
import           Data.Text.Lazy.Builder         ( Builder )
import           Lens.Micro

import           Mcmc.Trace

-- | Instruction about a parameter to monitor.
data MonitorBatch a = MonitorBatch
  {
    mpName :: Text               -- ^ Name of parameter.
  , mpFunc :: Trace a -> Builder -- ^ Instruction about how to extract the batch mean from the trace.
  }

-- TODO: Think about type of batch mean function. Probably Vector b -> b? ...
-- The trace type has to change!! But maybe, it is actually better to store the
-- parameters extracted by the monitors? I don't know.
--
-- It may be worth seeing if the implementation is too slow like this, and only
-- then, change.

-- | Batch monitor integral parameters such as 'Int'.
monitorBatchInt
  :: Integral b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchInt n l f = MonitorBatch n (T.decimal . f . map (^. l))

-- | Batch monitor real float parameters such as 'Double' with eight decimal
-- places (half precision).
monitorBatchRealFloat
  :: RealFloat b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchRealFloat n l f =
  MonitorBatch n (T.formatRealFloat T.Fixed (Just 8) . f . map (^. l))

-- | Batch monitor real float parameters such as 'Double' with full precision.
monitorBatchRealFloatF
  :: RealFloat b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchRealFloatF n l f = MonitorBatch n (T.realFloat . f . map (^. l))

-- | Batch monitor real float parameters such as 'Double' with scientific
-- notation and eight decimal places.
monitorBatchRealFloatS
  :: RealFloat b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchRealFloatS n l f =
  MonitorBatch n (T.formatRealFloat T.Exponent (Just 8) . f . map (^. l))
