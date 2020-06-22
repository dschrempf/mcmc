{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Monitor
-- Description :  Monitor a Markov chain
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu May 21 14:35:11 2020.
module Mcmc.Monitor
  ( -- * Create monitors
    Monitor (..),
    MonitorStdOut,
    monitorStdOut,
    MonitorFile,
    monitorFile,
    MonitorBatch,
    monitorBatch,

    -- * Use monitor
    mOpen,
    mHeader,
    mExec,
    mClose,
  )
where

import Data.Int
import qualified Data.Text.Lazy as T
import Data.Text.Lazy (Text)
import qualified Data.Text.Lazy.Builder as T
import qualified Data.Text.Lazy.IO as T
import Data.Time.Clock
import Mcmc.Item
import Mcmc.Monitor.Log
import Mcmc.Monitor.Parameter
import Mcmc.Monitor.ParameterBatch
import Mcmc.Monitor.Time
import Mcmc.Trace
import Numeric.Log
import System.IO
import Prelude hiding (sum)

-- | A 'Monitor' describes which part of the Markov chain should be logged and
-- where. Further, they allow output of summary statistics per iteration in a
-- flexible way.
data Monitor a
  = Monitor
      { -- | Monitor writing to standard output.
        mStdOut :: MonitorStdOut a,
        -- | Monitors writing to files.
        mFiles :: [MonitorFile a],
        -- | Monitors calculating batch means and
        -- writing to files.
        mBatches :: [MonitorBatch a]
      }

-- | Monitor to standard output.
data MonitorStdOut a
  = MonitorStdOut
      { msParams :: [MonitorParameter a],
        msPeriod :: Int
      }

-- | Monitor to standard output.
monitorStdOut ::
  -- | Instructions about which parameters to log.
  [MonitorParameter a] ->
  -- | Logging period.
  Int ->
  MonitorStdOut a
monitorStdOut ps p
  | p < 1 = error "monitorStdOut: Monitor period has to be 1 or larger."
  | otherwise = MonitorStdOut ps p

msIWidth :: Int64
msIWidth = 10

msWidth :: Int64
msWidth = 22

msRenderRow :: [Text] -> Text
msRenderRow xs = T.justifyRight msIWidth ' ' (head xs) <> T.concat vals
  where
    vals = map (T.justifyRight msWidth ' ') (tail xs)

msHeader :: MonitorStdOut a -> IO ()
msHeader m = T.hPutStr stdout $ T.unlines [row, sep]
  where
    row =
      msRenderRow $
        ["Iteration", "Log-Prior", "Log-Likelihood", "Log-Posterior"]
          ++ nms
          ++ ["Runtime", "ETA"]
    sep = T.replicate (T.length row) "─"
    nms = [T.pack $ mpName p | p <- msParams m]

msExec ::
  Int ->
  Item a ->
  UTCTime ->
  Int ->
  MonitorStdOut a ->
  IO ()
msExec i (Item x p l) st j m
  | i `mod` msPeriod m /= 0 =
    return ()
  | otherwise = do
      ct <- getCurrentTime
      let dt = ct `diffUTCTime` st
          timePerIter = dt / fromIntegral i
          eta = if i < 10
                then ""
                else renderDuration $ timePerIter * fromIntegral (j - i)
      T.hPutStrLn stdout
        $ msRenderRow
        $ [T.pack (show i), renderLog p, renderLog l, renderLog (p * l)]
          ++ [T.toLazyText $ mpFunc mp x | mp <- msParams m]
          ++ [renderDuration dt , eta]

-- | Monitor to a file.
data MonitorFile a
  = MonitorFile
      { mfName :: String,
        mfHandle :: Maybe Handle,
        mfParams :: [MonitorParameter a],
        mfPeriod :: Int
      }

-- XXX: The file monitor also includes iteration, prior, likelihood, and
-- posterior. What if I want to log trees; or other complex objects? In this
-- case, we need a simpler monitor to a file.

-- | Monitor parameters to a file.
monitorFile ::
  -- | Name; used as part of the file name.
  String ->
  -- | Instructions about which parameters to log.
  [MonitorParameter a] ->
  -- | Logging period.
  Int ->
  MonitorFile a
monitorFile n ps p
  | p < 1 = error "monitorFile: Monitor period has to be 1 or larger."
  | otherwise = MonitorFile n Nothing ps p

mfRenderRow :: [Text] -> Text
mfRenderRow = T.intercalate "\t"

mfOpen :: String -> MonitorFile a -> IO (MonitorFile a)
mfOpen n m = do
  h <- openFile (n <> mfName m <> ".monitor") WriteMode
  hSetBuffering h LineBuffering
  return $ m {mfHandle = Just h}

mfHeader :: MonitorFile a -> IO ()
mfHeader m = case mfHandle m of
  Nothing ->
    error $
      "mfHeader: No handle available for monitor with name "
        <> mfName m
        <> "."
  Just h ->
    T.hPutStrLn h
      $ mfRenderRow
      $ ["Iteration", "Log-Prior", "Log-Likelihood", "Log-Posterior"]
        ++ [T.pack $ mpName p | p <- mfParams m]

mfExec ::
  Int ->
  Item a ->
  MonitorFile a ->
  IO ()
mfExec i (Item x p l) m
  | i `mod` mfPeriod m /= 0 = return ()
  | otherwise = case mfHandle m of
    Nothing ->
      error $
        "mfExec: No handle available for monitor with name "
          <> mfName m
          <> "."
    Just h ->
      T.hPutStrLn h
        $ mfRenderRow
        $ T.pack (show i)
          : renderLog p
          : renderLog l
          : renderLog (p * l)
          : [T.toLazyText $ mpFunc mp x | mp <- mfParams m]

mfClose :: MonitorFile a -> IO ()
mfClose m = case mfHandle m of
  Just h -> hClose h
  Nothing -> error $ "mfClose: File was not opened for monitor " <> mfName m <> "."

-- | Monitor to a file, but calculate batch means for the given batch size.
--
-- XXX: Batch monitors are slow at the moment because the monitored parameter
-- has to be extracted from the state for each iteration.
data MonitorBatch a
  = MonitorBatch
      { mbName :: String,
        mbHandle :: Maybe Handle,
        mbParams :: [MonitorParameterBatch a],
        mbSize :: Int
      }

-- XXX: The batch monitor also includes iteration, prior, likelihood, and
-- posterior. What if I want to log trees; or other complex objects? In this
-- case, we need a simpler monitor to a file.

-- | Monitor parameters to a file, see 'MonitorBatch'.
monitorBatch ::
  -- | Name; used as part of the file name.
  String ->
  -- | Instructions about which parameters to log
  -- and how to calculate the batch means.
  [MonitorParameterBatch a] ->
  -- | Batch size.
  Int ->
  MonitorBatch a
monitorBatch n ps p
  | p < 2 = error "monitorBatch: Batch size has to be 2 or larger."
  | otherwise = MonitorBatch n Nothing ps p

mbOpen :: String -> MonitorBatch a -> IO (MonitorBatch a)
mbOpen n m = do
  h <- openFile (n <> mbName m <> ".batch") WriteMode
  hSetBuffering h LineBuffering
  return $ m {mbHandle = Just h}

mbHeader :: MonitorBatch a -> IO ()
mbHeader m = case mbHandle m of
  Nothing ->
    error $
      "mbHeader: No handle available for batch monitor with name "
        <> mbName m
        <> "."
  Just h ->
    T.hPutStrLn h
      $ mfRenderRow
      $ ["Iteration", "Mean log-Prior", "Mean log-Likelihood", "Mean log-Posterior"]
        ++ [T.pack $ mbpName mbp | mbp <- mbParams m]

logMean :: [Log Double] -> Log Double
logMean xs = sum xs / fromIntegral (length xs)

mbExec ::
  Int ->
  Trace a ->
  MonitorBatch a ->
  IO ()
mbExec i t' m
  | (i `mod` mbSize m /= 0) || (i == 0) = return ()
  | otherwise = case mbHandle m of
    Nothing ->
      error $
        "mbExec: No handle available for batch monitor with name "
          <> mbName m
          <> "."
    Just h ->
      T.hPutStrLn h
        $ mfRenderRow
        $ T.pack (show i)
          : renderLog mlps
          : renderLog mlls
          : renderLog mlos
          : [T.toLazyText $ mbpFunc mbp (map state t) | mbp <- mbParams m]
  where
    t = takeT (mbSize m) t'
    lps = map logPrior t
    lls = map logLikelihood t
    los = zipWith (*) lps lls
    mlps = logMean lps
    mlls = logMean lls
    mlos = logMean los

mbClose :: MonitorBatch a -> IO ()
mbClose m = case mbHandle m of
  Just h -> hClose h
  Nothing -> error $ "mfClose: File was not opened for batch monitor: " <> mbName m <> "."

-- | Open the files associated with the 'Monitor'.
mOpen :: String -> Monitor a -> IO (Monitor a)
mOpen n (Monitor s fs bs) = do
  fs' <- mapM (mfOpen n) fs
  mapM_ mfHeader fs'
  bs' <- mapM (mbOpen n) bs
  mapM_ mbHeader bs'
  return $ Monitor s fs' bs'

-- | Print header line of 'Monitor' (standard output only).
mHeader :: Monitor a -> IO ()
mHeader (Monitor s _ _) = msHeader s

-- | Execute monitors; print status information to standard output and files.
mExec ::
  -- | Iteration.
  Int ->
  -- | Start time.
  UTCTime ->
  -- | Trace of Markov chain.
  Trace a ->
  -- | Total number of iterations; to calculate ETA.
  Int ->
  -- | The monitor.
  Monitor a ->
  IO ()
mExec i t xs j (Monitor s fs bs) = do
  msExec i (headT xs) t j s
  mapM_ (mfExec i $ headT xs) fs
  mapM_ (mbExec i xs) bs

-- | Close the files associated with the 'Monitor'.
mClose :: Monitor a -> IO (Monitor a)
mClose m@(Monitor _ fms bms) = do
  mapM_ mfClose fms
  mapM_ mbClose bms
  let fms' = map (\fm -> fm {mfHandle = Nothing}) fms
  let bms' = map (\bm -> bm {mbHandle = Nothing}) bms
  return m {mFiles = fms', mBatches = bms'}
