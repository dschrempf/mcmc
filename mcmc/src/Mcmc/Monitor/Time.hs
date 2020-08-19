{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Monitor.Time
-- Description :  Print time related values
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Fri May 29 12:36:43 2020.
module Mcmc.Monitor.Time
  ( renderDuration,
    renderDurationS,
  )
where

import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Time.Clock
import Mcmc.Internal.ByteString

-- | Adapted from System.ProgressBar.renderDuration of package
-- [terminal-progressbar-0.4.1](https://hackage.haskell.org/package/terminal-progress-bar-0.4.1).
renderDuration :: NominalDiffTime -> BL.ByteString
renderDuration dt = hTxt <> mTxt <> sTxt
  where
    hTxt = renderDecimal h <> ":"
    mTxt = renderDecimal m <> ":"
    sTxt = renderDecimal s
    (h, hRem) = ts `quotRem` 3600
    (m, s) = hRem `quotRem` 60
    -- Total amount of seconds
    ts :: Int
    ts = round dt
    renderDecimal n = alignRightWith '0' 2 $ BB.toLazyByteString $ BB.intDec n

-- | Render duration in seconds.
renderDurationS :: NominalDiffTime -> BL.ByteString
renderDurationS dt = BB.toLazyByteString $ BB.intDec ts
  where
    ts :: Int
    ts = round dt
