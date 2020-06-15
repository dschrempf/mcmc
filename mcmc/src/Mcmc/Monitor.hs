{-# LANGUAGE OverloadedStrings #-}

{- |
Module      :  Mcmc.Monitor
Description :  Monitor a Markov chain
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 21 14:35:11 2020.

-}

-- TODO Monitor every iteration by default and note that subsampling is
-- disadvantageous with reference to Geyer.

-- TODO Allow batch mean monitors. Allow custom function to calculate the means.

module Mcmc.Monitor
  (
    -- * Create monitors
    Monitor(..)
  , MonitorStdOut
  , monitorStdOut
  , MonitorFile
  , monitorFile
    -- * Use monitor
  , mOpen
  , mHeader
  , mExec
  , mClose
  )
where

import           Data.Int
import qualified Data.Text.Lazy                as T
import qualified Data.Text.Lazy.IO             as T
import qualified Data.Text.Lazy.Builder        as T
import           Data.Text.Lazy                 ( Text )
import           Data.Time.Clock
import           Numeric.Log
import           System.IO

import           Mcmc.Monitor.Log
import           Mcmc.Monitor.Parameter
import           Mcmc.Monitor.Time

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
monitorStdOut ps p | p < 1    = error "monitorStdOut: Monitor period has to be 1 or larger."
                   | otherwise = MonitorStdOut ps p

msIWidth :: Int64
msIWidth = 10

msWidth :: Int64
msWidth = 22

msRenderRow :: [Text] -> Text
msRenderRow xs = T.justifyRight msIWidth ' ' (head xs) <> T.concat vals
  where vals = map (T.justifyRight msWidth ' ') (tail xs)

msHeader :: MonitorStdOut a -> IO ()
msHeader m = T.hPutStr stdout $ T.unlines [row, sep]
 where
  row =
    msRenderRow
      $  ["Iteration", "Log-Prior", "Log-Likelihood", "Log-Posterior"]
      ++ nms
      ++ ["Runtime", "ETA"]
  sep = T.replicate (T.length row) "─"
  nms = [ mpName p | p <- msParams m ]

msExec
  :: Int
  -> Log Double
  -> Log Double
  -> Log Double
  -> NominalDiffTime
  -> a
  -> Int
  -> MonitorStdOut a
  -> IO ()
msExec i p l o t x j m
  | i `mod` msPeriod m /= 0
  = return ()
  | otherwise
  = T.hPutStrLn stdout
    $  msRenderRow
    $  [T.pack (show i), renderLog p, renderLog l, renderLog o]
    ++ [ T.toLazyText $ mpFunc mp x | mp <- msParams m ]
    ++ [renderDuration t, eta]
 where
  eta =
    if i < 10 then "" else renderDuration $ timePerIter * fromIntegral (j - i)
  timePerIter = t / fromIntegral i

-- | Monitor to a file.
data MonitorFile a = MonitorFile
  {
    mfFile   :: FilePath
  , mfHandle :: Maybe Handle
  , mfParams :: [MonitorParameter a]
  , mfPeriod :: Int
  }

-- TODO: This monitor also includes iteration, prior, likelihood, and posterior.
-- What if I want to log trees; or other complex objects? In this case, we need
-- a simpler monitor to a file.

-- | Monitor writing to a file.
monitorFile
  :: FilePath             -- ^ File path; file will be overwritten!
  -> [MonitorParameter a] -- ^ Instructions about which parameters to log.
  -> Int                  -- ^ Logging period.
  -> MonitorFile a
monitorFile f ps p | p <= 1    = error "monitorFile: Monitor period has to be 1 or larger."
                   | otherwise = MonitorFile f Nothing ps p

mfRenderRow :: [Text] -> Text
mfRenderRow = T.intercalate "\t"

mfOpen :: MonitorFile a -> IO (MonitorFile a)
mfOpen m = do
  h <- openFile (mfFile m) WriteMode
  return $ m { mfHandle = Just h }

mfHeader :: MonitorFile a -> IO ()
mfHeader m = case mfHandle m of
  Nothing ->
    error
      $  "mfHeader: No handle available for monitor with file "
      <> mfFile m
      <> "."
  Just h ->
    T.hPutStrLn h
      $  mfRenderRow
      $  ["Iteration", "Log-Prior", "Log-Likelihood", "Log-Posterior"]
      ++ [ mpName p | p <- mfParams m ]

mfExec
  :: Int
  -> Log Double
  -> Log Double
  -> Log Double
  -> a
  -> MonitorFile a
  -> IO ()
mfExec i p l o x m
  | i `mod` mfPeriod m /= 0 = return ()
  | otherwise = case mfHandle m of
    Nothing ->
      error
        $  "mfExec: No handle available for monitor with file "
        <> mfFile m
        <> "."
    Just h ->
      T.hPutStrLn h
        $ mfRenderRow
        $ T.pack (show i)
        : renderLog p
        : renderLog l
        : renderLog o
        : [ T.toLazyText $ mpFunc mp x | mp <- mfParams m ]

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

-- | Execute monitors; print status information to standard output and files.
mExec
  :: Int -- ^ Iteration.
  -> Log Double -- ^ Prior.
  -> Log Double -- ^ Likelihood.
  -> NominalDiffTime -- ^ Run time.
  -> a -- ^ State of Markov chain.
  -> Int -- ^ Total number of iterations; to calculate ETA.
  -> Monitor a -- ^ The monitor.
  -> IO ()
mExec i p l t x j (Monitor s fs) =
  msExec i p l (p * l) t x j s >> mapM_ (mfExec i p l (p * l) x) fs

-- | Close the files associated with the 'Monitor'.
mClose :: Monitor a -> IO (Monitor a)
mClose m@(Monitor _ fms) = do
  mapM_ mfClose fms
  let fms' = map (\fm -> fm { mfHandle = Nothing }) fms
  return m { mFiles = fms' }

