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
  ( MonitorBatch(..)
  , monitorBatchInt
  , monitorBatchRealFloat
  , monitorBatchRealFloatF
  , monitorBatchRealFloatS
  )
where

import qualified Data.Text.Lazy.Builder.Int    as T
import qualified Data.Text.Lazy.Builder.RealFloat
                                               as T
import           Data.Text.Lazy                 ( Text )
import           Data.Text.Lazy.Builder         ( Builder )
import           Lens.Micro

import           Mcmc.Trace

-- | Instruction about a parameter to monitor via batch means. Usually, the
-- monitored parameter is average over the batch size. However, arbitrary
-- functions performing more complicated analyses on the states in the batch can
-- be provided.
--
-- XXX: Slow at the moment because the monitored parameter has to be extracted
-- from the state for each iteration.
--
-- XXX: Even slower at the moment because normal lists are used.
data MonitorBatch a = MonitorBatch
  {
    mpName :: Text               -- ^ Name of parameter.
  , mpFunc :: Trace a -> Builder -- ^ Instruction about how to extract the batch mean from the trace.
  }

params :: Lens' a b -> [a] -> [b]
params l = map (^. l)

-- | Batch monitor integral parameters such as 'Int'.
monitorBatchInt
  :: Integral b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchInt n l f = MonitorBatch n (T.decimal . f . params l . states)

-- | Batch monitor real float parameters such as 'Double' with eight decimal
-- places (half precision).
monitorBatchRealFloat
  :: RealFloat b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchRealFloat n l f =
  MonitorBatch n (T.formatRealFloat T.Fixed (Just 8) . f . params l . states)

-- | Batch monitor real float parameters such as 'Double' with full precision.
monitorBatchRealFloatF
  :: RealFloat b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchRealFloatF n l f = MonitorBatch n (T.realFloat . f . params l . states)

-- | Batch monitor real float parameters such as 'Double' with scientific
-- notation and eight decimal places.
monitorBatchRealFloatS
  :: RealFloat b
  => Text       -- ^ Name of monitor.
  -> Lens' a b  -- ^ Instruction about which parameter to monitor.
  -> ([b] -> b) -- ^ Function to calculate the batch mean.
  -> MonitorBatch a
monitorBatchRealFloatS n l f =
  MonitorBatch n (T.formatRealFloat T.Exponent (Just 8) . f . params l . states)
