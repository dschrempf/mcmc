{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Monitor.Parameter
-- Description :  Monitor parameters
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 11:11:49 2020.
module Mcmc.Monitor.Parameter
  ( -- * Parameter monitors
    MonitorParameter (..),
    (>$<),
    monitorInt,
    monitorDouble,
    monitorDoubleF,
    monitorDoubleE,
  )
where

import qualified Data.ByteString.Builder as BB
import qualified Data.Double.Conversion.ByteString as BC
import Data.Functor.Contravariant

-- | Instruction about a parameter to monitor.
--
-- Convert a parameter monitor from one data type to another with '(>$<)'.
--
-- For example, monitor a 'Double' value being the first entry of a tuple:
--
-- @
-- mon = fst >$< monitorDouble
-- @
data MonitorParameter a = MonitorParameter
  { -- | Name of parameter.
    mpName :: String,
    -- | Instruction about how to extract the parameter from the state.
    mpFunc :: a -> BB.Builder
  }

instance Contravariant MonitorParameter where
  contramap f (MonitorParameter n m) = MonitorParameter n (m . f)

-- | Monitor 'Int'.
monitorInt ::
  -- | Name.
  String ->
  MonitorParameter Int
monitorInt n = MonitorParameter n BB.intDec

-- | Monitor 'Double' with eight decimal places (half precision).
monitorDouble ::
  -- | Name.
  String ->
  MonitorParameter Double
monitorDouble n = MonitorParameter n (BB.byteString . BC.toFixed 8)

-- | Monitor 'Double' with full precision computing the shortest string of
-- digits that correctly represent the number.
monitorDoubleF ::
  -- | Name.
  String ->
  MonitorParameter Double
monitorDoubleF n = MonitorParameter n (BB.byteString . BC.toShortest)

-- | Monitor 'Double' in exponential format and half precision.
monitorDoubleE ::
  -- | Name.
  String ->
  MonitorParameter Double
monitorDoubleE n = MonitorParameter n (BB.byteString . BC.toExponential 8)
