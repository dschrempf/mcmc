{-# LANGUAGE OverloadedStrings #-}

{- |
Module      :  Mcmc.Monitor.Time
Description :  Print time related values
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Fri May 29 12:36:43 2020.

-}

module Mcmc.Monitor.Time
  ( renderDuration
  ) where

import qualified Data.Text.Lazy             as T
import qualified Data.Text.Lazy.Builder     as T
import qualified Data.Text.Lazy.Builder.Int as T
import Data.Text.Lazy (Text)
import Data.Time.Clock

-- | Adapted from System.ProgressBar.renderDuration of package
-- [terminal-progressbar-0.4.1](https://hackage.haskell.org/package/terminal-progress-bar-0.4.1).
renderDuration :: NominalDiffTime -> Text
renderDuration dt = hTxt <> mTxt <> sTxt
  where
    hTxt = renderDecimal h <> ":"
    mTxt = renderDecimal m <> ":"
    sTxt = renderDecimal s

    (h, hRem) = ts   `quotRem` 3600
    (m, s   ) = hRem `quotRem`   60

    -- Total amount of seconds
    ts :: Int
    ts = round dt

    renderDecimal n = T.justifyRight 2 '0' $ T.toLazyText $ T.decimal n
