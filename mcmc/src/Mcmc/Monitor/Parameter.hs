{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Mcmc.Monitor.Parameter
Description :  Monitor parameters
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May 29 11:11:49 2020.

-}

module Mcmc.Monitor.Parameter
  ( MonitorParameter(..)
  , monitorInt
  , monitorRealFloat
  , monitorRealFloatF
  , monitorRealFloatS
  )
where

-- import qualified Data.Text.Lazy         as T
import qualified Data.Text.Lazy.Builder.Int    as T
import qualified Data.Text.Lazy.Builder.RealFloat
                                               as T
import           Data.Text.Lazy                 ( Text )
import           Data.Text.Lazy.Builder         ( Builder )
import           Lens.Micro

-- | Instruction about a parameter to monitor.
data MonitorParameter a = MonitorParameter
  {
    mpName :: Text         -- ^ Name of parameter.
  , mpFunc :: a -> Builder -- ^ Instruction about how to extract the parameter
                           -- from the state.
  }

-- | Monitor integral parameters such as 'Int'.
monitorInt
  :: Integral b
  => Text -- ^ Name of monitor.
  -> Lens' a b
  -> MonitorParameter a
monitorInt n l = MonitorParameter n (\x -> T.decimal $ x ^. l)

-- | Monitor real float parameters such as 'Double' with eight decimal places
-- (half precision).
monitorRealFloat
  :: RealFloat b
  => Text -- ^ Name of monitor.
  -> Lens' a b
  -> MonitorParameter a
monitorRealFloat n l =
  MonitorParameter n (\x -> T.formatRealFloat T.Fixed (Just 8) $ x ^. l)

-- | Monitor real float parameters such as 'Double' with full precision.
monitorRealFloatF
  :: RealFloat b
  => Text -- ^ Name of monitor.
  -> Lens' a b
  -> MonitorParameter a
monitorRealFloatF n l = MonitorParameter n (\x -> T.realFloat $ x ^. l)

-- | Monitor real float parameters such as 'Double' with scientific notation and
-- eight decimal places.
monitorRealFloatS
  :: RealFloat b
  => Text -- ^ Name of monitor.
  -> Lens' a b
  -> MonitorParameter a
monitorRealFloatS n l =
  MonitorParameter n (\x -> T.formatRealFloat T.Exponent (Just 8) $ x ^. l)
