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
  ( MonitorParameter (..),
    monitorInt,
    monitorRealFloat,
    monitorRealFloatF,
    monitorRealFloatS,
  )
where

import Data.Text.Lazy (Text)
import Data.Text.Lazy.Builder (Builder)
import qualified Data.Text.Lazy.Builder.Int as T
import qualified Data.Text.Lazy.Builder.RealFloat as T
import Lens.Micro

-- | Instruction about a parameter to monitor.
data MonitorParameter a
  = MonitorParameter
      { -- | Name of parameter.
        mpName :: Text,
        -- | Instruction about how to extract the parameter from the state.
        mpFunc :: a -> Builder
      }

-- | Monitor integral parameters such as 'Int'.
monitorInt ::
  Integral b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameter a
monitorInt n l = MonitorParameter n (\x -> T.decimal $ x ^. l)
{-# SPECIALIZE monitorInt :: Text -> Lens' a Int -> MonitorParameter a #-}

-- | Monitor real float parameters such as 'Double' with eight decimal places
-- (half precision).
monitorRealFloat ::
  RealFloat b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameter a
monitorRealFloat n l =
  MonitorParameter n (\x -> T.formatRealFloat T.Fixed (Just 8) $ x ^. l)
{-# SPECIALIZE monitorRealFloat :: Text -> Lens' a Double -> MonitorParameter a #-}

-- | Monitor real float parameters such as 'Double' with full precision.
monitorRealFloatF ::
  RealFloat b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameter a
monitorRealFloatF n l = MonitorParameter n (\x -> T.realFloat $ x ^. l)
{-# SPECIALIZE monitorRealFloatF :: Text -> Lens' a Double -> MonitorParameter a #-}

-- | Monitor real float parameters such as 'Double' with scientific notation and
-- eight decimal places.
monitorRealFloatS ::
  RealFloat b =>
  -- | Name of monitor.
  Text ->
  -- | Instruction about which parameter to monitor.
  Lens' a b ->
  MonitorParameter a
monitorRealFloatS n l =
  MonitorParameter
    n
    (\x -> T.formatRealFloat T.Exponent (Just 8) $ x ^. l)
{-# SPECIALIZE monitorRealFloatS :: Text -> Lens' a Double -> MonitorParameter a #-}
