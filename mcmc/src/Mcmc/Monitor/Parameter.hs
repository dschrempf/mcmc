{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc.Monitor.Parameter
-- Description :  Monitor parameters
-- Copyright   :  (c) Dominik Schrempf, 2020
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
    (@.),
    monitorInt,
    monitorRealFloat,
    monitorRealFloatF,
    monitorRealFloatS,
  )
where

import Data.Text.Lazy.Builder (Builder)
import qualified Data.Text.Lazy.Builder.Int as T
import qualified Data.Text.Lazy.Builder.RealFloat as T
import Lens.Micro

-- | Instruction about a parameter to monitor.
data MonitorParameter a
  = MonitorParameter
      { -- | Name of parameter.
        mpName :: String,
        -- | Instruction about how to extract the parameter from the state.
        mpFunc :: a -> Builder
      }

-- | Convert a parameter monitor from one data type to another using a lens.
--
-- For example, to monitor a real float value being the first entry of a tuple:
--
-- @
-- mon = _1 @@ monitorRealFloat
-- @
(@.) :: Lens' b a -> MonitorParameter a -> MonitorParameter b
(@.) l (MonitorParameter n f) = MonitorParameter n (\x -> f $ x^.l)

-- | Monitor integral parameters such as 'Int'.
monitorInt ::
  Integral a =>
  -- | Name of monitor.
  String ->
  MonitorParameter a
monitorInt n = MonitorParameter n T.decimal
{-# SPECIALIZE monitorInt :: String -> MonitorParameter Int #-}

-- | Monitor real float parameters such as 'Double' with eight decimal places
-- (half precision).
monitorRealFloat ::
  RealFloat a =>
  -- | Name of monitor.
  String ->
  MonitorParameter a
monitorRealFloat n = MonitorParameter n (T.formatRealFloat T.Fixed (Just 8))
{-# SPECIALIZE monitorRealFloat :: String -> MonitorParameter Double #-}

-- | Monitor real float parameters such as 'Double' with full precision.
monitorRealFloatF ::
  RealFloat a =>
  -- | Name of monitor.
  String ->
  MonitorParameter a
monitorRealFloatF n = MonitorParameter n T.realFloat
{-# SPECIALIZE monitorRealFloatF :: String -> MonitorParameter Double #-}

-- | Monitor real float parameters such as 'Double' with scientific notation and
-- eight decimal places.
monitorRealFloatS ::
  RealFloat a =>
  -- | Name of monitor.
  String ->
  MonitorParameter a
monitorRealFloatS n = MonitorParameter n (T.formatRealFloat T.Exponent (Just 8))
{-# SPECIALIZE monitorRealFloatS :: String -> MonitorParameter Double #-}
