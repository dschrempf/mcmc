-- |
-- Module      :  Mcmc.Internal.ByteString
-- Description :  ByteString tools
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Aug  3 10:46:27 2020.
module Mcmc.Internal.ByteString
  ( alignRightWith,
    alignRight,
    alignLeftWith,
    alignLeft,
  )
where

import qualified Data.ByteString.Lazy.Char8 as BL

-- | For a given width, align string to the right; use given fill character;
-- trim on the left if string is longer.
alignRightWith :: Char -> Int -> BL.ByteString -> BL.ByteString
alignRightWith c n s =
  BL.replicate (fromIntegral n - l) c <> BL.take (fromIntegral n) s
  where
    l = BL.length s

-- | For a given width, align string to the right; trim on the left if string is
-- longer.
alignRight :: Int -> BL.ByteString -> BL.ByteString
alignRight = alignRightWith ' '

-- | For a given width, align string to the left; use given fill character; trim
-- on the right if string is longer.
alignLeftWith :: Char -> Int -> BL.ByteString -> BL.ByteString
alignLeftWith c n s =
  BL.take (fromIntegral n) s <> BL.replicate (fromIntegral n - l) c
  where
    l = BL.length s

-- | For a given width, align string to the left; trim on the right if string is
-- longer.
alignLeft :: Int -> BL.ByteString -> BL.ByteString
alignLeft = alignLeftWith ' '
