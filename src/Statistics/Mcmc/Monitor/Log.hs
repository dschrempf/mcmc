{- |
Module      :  Statistics.Mcmc.Monitor.Log
Description :  Monitor logarithmic values
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May 29 12:30:03 2020.

-}

module Statistics.Mcmc.Monitor.Log
  ( renderLog
  ) where

import qualified Data.Text.Lazy.Builder           as T
import qualified Data.Text.Lazy.Builder.RealFloat as T
import Data.Text.Lazy (Text)
import Numeric.Log

-- | Print a log value.
renderLog :: Log Double -> Text
renderLog = T.toLazyText . T.formatRealFloat T.Fixed (Just 8) . ln
{-# INLINEABLE renderLog  #-}
