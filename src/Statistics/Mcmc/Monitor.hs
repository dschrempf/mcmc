{-# LANGUAGE OverloadedStrings #-}

{- |
Module      :  Statistics.Mcmc.Monitor
Description :  Monitor a Markov chain
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 21 14:35:11 2020.

-}

module Statistics.Mcmc.Monitor
  (
    -- * Create monitors
    Monitor (..)
  , MonitorStdOut
  , monitorStdOut
  , MonitorFile
  , monitorFile
    -- * Use monitor
  , mOpen
  , mHeader
  , mExec
  , mClose
  ) where

-- TODO: Use builders (easier for users).

import Data.Int
import qualified Data.List.NonEmpty     as N
import Data.List.NonEmpty (NonEmpty)
import qualified Data.Text.Lazy         as T
import qualified Data.Text.Lazy.IO      as T
import qualified Data.Text.Lazy.Builder as T
import Data.Text.Lazy (Text)
import System.IO

import Statistics.Mcmc.Monitor.Parameter

-- | A 'Monitor' describes which part of the Markov chain should be logged and
-- where. Further, they allow output of summary statistics per iteration in a
-- flexible way.
data Monitor a = Monitor
 {
   mStdOut :: MonitorStdOut a -- ^ Monitor writing to standard output.
 , mFiles  :: [MonitorFile a] -- ^ Monitors writing to files.
 }

-- | Monitor to standard output.
data MonitorStdOut a = MonitorStdOut
  {
    msParams :: [MonitorParameter a]
  , msPeriod :: Int
  }

-- | Monitor to standard output.
monitorStdOut
  :: [MonitorParameter a] -- ^ Instructions about which parameters to log.
  -> Int                  -- ^ Logging period.
  -> MonitorStdOut a
monitorStdOut = MonitorStdOut

msIWidth :: Int64
msIWidth = 14

msWidth :: Int64
msWidth = 22

msRenderRow :: NonEmpty Text -> Text
msRenderRow xs = T.justifyRight msIWidth ' ' (N.head xs) <> T.concat vals
  where
    vals = map (T.justifyRight msWidth ' ') (N.tail xs)

msHeader :: MonitorStdOut a -> IO ()
msHeader m = T.hPutStr stdout $ T.unlines [row, sep]
  where row = msRenderRow $ T.pack "Iteration" N.:| [ mpName p | p <- msParams m ]
        sep = T.replicate (T.length row) "─"

msExec :: Int -> a -> MonitorStdOut a -> IO ()
msExec i x m
  | i `mod` msPeriod m /= 0        = return ()
  | otherwise                    = row
  where row = T.hPutStrLn stdout $ msRenderRow $
          T.pack (show i) N.:| [ T.toLazyText $ mpFunc p x | p <- msParams m ]

-- | Monitor to a file.
data MonitorFile a = MonitorFile
  {
    mfFile   :: FilePath
  , mfHandle :: Maybe Handle
  , mfParams :: [MonitorParameter a]
  , mfPeriod :: Int
  }

-- | Monitor writing to a file.
monitorFile
  :: FilePath             -- ^ File path; file will be overwritten!
  -> [MonitorParameter a] -- ^ Instructions about which parameters to log.
  -> Int                  -- ^ Logging period.
  -> MonitorFile a
monitorFile f = MonitorFile f Nothing

mfRenderRow :: [Text] -> Text
mfRenderRow = T.intercalate "\t"

mfOpen :: MonitorFile a -> IO (MonitorFile a)
mfOpen m = do
  h <- openFile (mfFile m) WriteMode
  return $ m { mfHandle = Just h }

mfHeader :: MonitorFile a -> IO ()
mfHeader m =
  case mfHandle m of
    Nothing -> error $ "mfHeader: No handle available for monitor with file " <> mfFile m <> "."
    Just h  -> T.hPutStrLn h row
  where
    row = mfRenderRow $ T.pack "Iteration" : [ mpName p | p <- mfParams m ]

mfExec :: Int -> a -> MonitorFile a -> IO ()
mfExec i x m
  | i `mod` mfPeriod m /= 0 = return ()
  | otherwise =
    case mfHandle m of
      Nothing -> error $ "mfExec: No handle available for monitor with file " <> mfFile m <> "."
      Just h  -> T.hPutStrLn h $ mfRenderRow $
        T.pack (show i) : [T.toLazyText $ mpFunc p x | p <- mfParams m ]

mfClose :: MonitorFile a -> IO ()
mfClose m = case mfHandle m of
  Just h  -> hClose h
  Nothing -> error $ "mfClose: File was not opened: " <> mfFile m <> "."

-- | Open the files associated with the 'Monitor'.
mOpen :: Monitor a -> IO (Monitor a)
mOpen (Monitor s fs) = do
  fs' <- mapM mfOpen fs
  mapM_ mfHeader fs'
  return $ Monitor s fs'

-- | Print header line of 'Monitor' (standard output only).
mHeader :: Monitor a -> IO ()
mHeader (Monitor s _) = msHeader s

-- TODO: Provide and monitor prior, likelihood, posterior.

-- TODO: Provide and monitor run time and ETA.

-- | Print logs for given iteration.
mExec :: Int -> a -> Monitor a -> IO ()
mExec i x (Monitor s fs) = msExec i x s >> mapM_ (mfExec i x) fs

-- | Close the files associated with the 'Monitor'.
mClose :: Monitor a -> IO ()
mClose (Monitor _ fs) = mapM_ mfClose fs
