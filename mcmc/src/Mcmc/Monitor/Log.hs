-- |
-- Module      :  Mcmc.Monitor.Log
-- Description :  Monitor logarithmic values
-- Copyright   :  (c) Dominik Schrempf, 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 12:30:03 2020.
module Mcmc.Monitor.Log
  ( renderLog,
  )
where

import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Double.Conversion.ByteString as BC
import Numeric.Log

-- | Print a log value.
renderLog :: Log Double -> BL.ByteString
renderLog = BL.fromStrict . BC.toFixed 8 . ln
{-# INLINEABLE renderLog #-}
