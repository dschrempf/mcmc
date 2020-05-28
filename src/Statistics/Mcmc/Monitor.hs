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

-- TODO: Separate burn in with auto tune, list of moves, and monitors.

module Statistics.Mcmc.Monitor
  (
    -- * Create monitors
    Monitor (..)
  , MonitorParameter (..)
  , MonitorSimple (..)
  , MonitorFile
  , monitorFile
    -- * Use monitor
  , mOpen
  , mHeader
  , mExec
  , mClose
  ) where

import qualified Data.List.NonEmpty as N
import Data.List.NonEmpty (NonEmpty)
import qualified Data.Text          as T
import qualified Data.Text.IO       as T
import Data.Text (Text)
import System.IO

-- | A 'Monitor' describes which part of the Markov chain should be logged and
-- where. Further, they allow output of summary statistics per iteration in a
-- flexible way.
data Monitor a = Monitor
 {
   mStdOut :: MonitorSimple a -- ^ Monitor writing to standard output.
 , mFiles  :: [MonitorFile a] -- ^ Monitors writing to files.
 }

-- | Instruction about a parameter to monitor.
data MonitorParameter a = MonitorParameter
  {
    mpName :: Text      -- ^ Name of parameter.
  , mpFunc :: a -> Text -- ^ Instruction about how to extract the parameter from
                        -- the state.
  }

-- | Monitor a variable of the state space. The key part is the function 'mShow'
-- which describes what should be logged. Several monitors for standard types
-- are provided.
data MonitorSimple a = MonitorSimple
  {
    msParams :: [MonitorParameter a] -- ^ Parameters to monitor.
  , msFreq   :: Int                  -- ^ Logging period.
  }

renderRow :: NonEmpty Text -> Text
renderRow xs = T.justifyRight 6 ' ' (N.head xs) <> T.concat vals
  where vals = map (T.justifyRight 10 ' ') (N.tail xs)

msHeader :: MonitorSimple a -> Handle -> IO ()
msHeader m h = T.hPutStr h $ T.unlines [row, T.replicate (T.length row) "─"]
  where row = renderRow $ T.pack "Cycle" N.:| [ mpName p | p <- msParams m ]

msExec :: Int -> a -> MonitorSimple a -> Handle -> IO ()
msExec i x m h
  | i `mod` msFreq m /= 0 = return ()
  | otherwise = T.hPutStrLn h $ renderRow $ T.pack (show i) N.:| [ mpFunc p x | p <- msParams m ]

-- | Monitor to a file.
data MonitorFile a = MonitorFile
  {
    mfFile   :: FilePath
  , mfHandle :: Maybe Handle
  , mfSimple :: MonitorSimple a
  }

-- | Monitor writing to a file.
monitorFile
  :: FilePath          -- ^ File path; file will be overwritten!
  -> [MonitorParameter a]       -- ^ Instructions about which parameters to log.
  -> Int               -- ^ Logging period.
  -> MonitorFile a
monitorFile f bs p = MonitorFile f Nothing (MonitorSimple bs p)

mfOpen :: MonitorFile a -> IO (MonitorFile a)
mfOpen m = do
  h <- openFile (mfFile m) WriteMode
  return $ m { mfHandle = Just h }

mfHeader :: MonitorFile a -> IO ()
mfHeader m = case mfHandle m of
                 Nothing -> error $ "mfExec: No handle available for monitor with file " <> mfFile m <> "."
                 Just h  -> msHeader (mfSimple m) h

mfExec :: Int -> a -> MonitorFile a -> IO ()
mfExec i x m = case mfHandle m of
                 Nothing -> error $ "mfExec: No handle available for monitor with file " <> mfFile m <> "."
                 Just h  -> msExec i x (mfSimple m) h

mfClose :: MonitorFile a -> IO ()
mfClose m = case mfHandle m of
  Just h  -> hClose h
  Nothing -> error $ "mfClose: File was not opened: " <> mfFile m <> "."

-- | Open the files associated with the 'Monitor'.
mOpen :: Monitor a -> IO (Monitor a)
mOpen (Monitor s fs) = do
  fs' <- mapM mfOpen fs
  return $ Monitor s fs'

-- | Print header lines of 'Monitor'.
mHeader :: Monitor a -> IO ()
mHeader (Monitor s fs) = msHeader s stdout >> mapM_ mfHeader fs

-- | Print logs for given iteration.
mExec :: Int -> a -> Monitor a -> IO ()
mExec i x (Monitor s fs) = msExec i x s stdout >> mapM_ (mfExec i x) fs

-- | Close the files associated with the 'Monitor'.
mClose :: Monitor a -> IO ()
mClose (Monitor _ fs) = mapM_ mfClose fs
