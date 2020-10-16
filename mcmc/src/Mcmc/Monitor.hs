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

    -- * Use monitors
    mOpen,
    mAppend,
    mExec,
    mClose,
  )
where

import Control.Monad
import qualified Data.ByteString.Builder as BB
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.Int
import Data.Time.Clock
import Mcmc.Internal.ByteString
import Mcmc.Item
import Mcmc.Monitor.Log
import Mcmc.Monitor.Parameter
import Mcmc.Monitor.ParameterBatch
import Mcmc.Monitor.Time
import Mcmc.Trace
import Mcmc.Verbosity
import Numeric.Log
import System.Directory
import System.IO
import Prelude hiding (sum)

-- | A 'Monitor' describes which part of the Markov chain should be logged and
-- where. Further, they allow output of summary statistics per iteration in a
-- flexible way.
data Monitor a = Monitor
  { -- | Monitor writing to standard output.
    mStdOut :: MonitorStdOut a,
    -- | Monitors writing to files.
    mFiles :: [MonitorFile a],
    -- | Monitors calculating batch means and
    -- writing to files.
    mBatches :: [MonitorBatch a]
  }

-- | Monitor to standard output; constructed with 'monitorStdOut'.
data MonitorStdOut a = MonitorStdOut
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

msIWidth :: Int
msIWidth = 12

msWidth :: Int
msWidth = 22

msRenderRow :: [BL.ByteString] -> BL.ByteString
msRenderRow xs = alignRight msIWidth (head xs) <> BL.concat vals
  where
    vals = map (alignRight msWidth) (tail xs)

msHeader :: MonitorStdOut a -> BL.ByteString
msHeader m = BL.intercalate "\n" [row, sep]
  where
    row =
      msRenderRow $
        ["Iteration", "Log-Prior", "Log-Likelihood", "Log-Posterior"]
          ++ nms
          ++ ["Runtime", "ETA"]
    sep = "   " <> BL.replicate (BL.length row - 3) '-'
    nms = [BL.pack $ mpName p | p <- msParams m]

msDataLine ::
  Int ->
  Item a ->
  Int ->
  UTCTime ->
  Int ->
  MonitorStdOut a ->
  IO BL.ByteString
msDataLine i (Item x p l) ss st j m = do
  ct <- getCurrentTime
  let dt = ct `diffUTCTime` st
      -- Careful, don't evaluate this when i == ss.
      timePerIter = dt / fromIntegral (i - ss)
      -- -- Always 0; doesn't make much sense.
      -- tpi = if (i - ss) < 10
      --       then ""
      --       else renderDurationS timePerIter
      eta =
        if (i - ss) < 10
          then ""
          else renderDuration $ timePerIter * fromIntegral (j - i)
  return $
    msRenderRow $
      [BL.pack (show i), renderLog p, renderLog l, renderLog (p * l)]
        ++ [BB.toLazyByteString $ mpFunc mp x | mp <- msParams m]
        ++ [renderDuration dt, eta]

msExec ::
  Int ->
  Item a ->
  Int ->
  UTCTime ->
  Int ->
  MonitorStdOut a ->
  IO (Maybe BL.ByteString)
msExec i it ss st j m
  | i `mod` msPeriod m /= 0 = return Nothing
  | i `mod` (msPeriod m * 100) == 0 = do
    l <- msDataLine i it ss st j m
    return $ Just $ msHeader m <> "\n" <> l
  | otherwise = Just <$> msDataLine i it ss st j m

-- | Monitor to a file; constructed with 'monitorFile'.
data MonitorFile a = MonitorFile
  { mfName :: String,
    mfHandle :: Maybe Handle,
    mfParams :: [MonitorParameter a],
    mfPeriod :: Int
  }

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

mfRenderRow :: [BL.ByteString] -> BL.ByteString
mfRenderRow = BL.intercalate "\t"

open' :: String -> Bool -> IO Handle
open' n frc = do
  fe <- doesFileExist n
  case (fe, frc) of
    (False, _) -> openFile n WriteMode
    (True, True) -> openFile n WriteMode
    (True, False) -> error $ "open': File \"" <> n <> "\" exists; probably use 'force'?"

mfOpen :: String -> Bool -> MonitorFile a -> IO (MonitorFile a)
mfOpen n frc m = do
  let mfn = n <> mfName m <> ".monitor"
  h <- open' mfn frc
  hSetBuffering h LineBuffering
  return $ m {mfHandle = Just h}

mfAppend :: String -> MonitorFile a -> IO (MonitorFile a)
mfAppend n m = do
  let fn = n <> mfName m <> ".monitor"
  fe <- doesFileExist fn
  if fe
    then do
      h <- openFile fn AppendMode
      hSetBuffering h LineBuffering
      return $ m {mfHandle = Just h}
    else error $ "mfAppend: Monitor file does not exist: " ++ fn ++ "."

mfHeader :: MonitorFile a -> IO ()
mfHeader m = case mfHandle m of
  Nothing ->
    error $
      "mfHeader: No handle available for monitor with name "
        <> mfName m
        <> "."
  Just h ->
    BL.hPutStrLn h $
      mfRenderRow $
        ["Iteration", "Log-Prior", "Log-Likelihood", "Log-Posterior"]
          ++ [BL.pack $ mpName p | p <- mfParams m]

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
      BL.hPutStrLn h $
        mfRenderRow $
          BL.pack (show i) :
          renderLog p :
          renderLog l :
          renderLog (p * l) :
            [BB.toLazyByteString $ mpFunc mp x | mp <- mfParams m]

mfClose :: MonitorFile a -> IO ()
mfClose m = case mfHandle m of
  Just h -> hClose h
  Nothing -> error $ "mfClose: File was not opened for monitor " <> mfName m <> "."

-- | Monitor to a file, but calculate batch means for the given batch size;
-- constructed with 'monitorBatch'.
--
-- Batch monitors are slow at the moment because the monitored parameter has to
-- be extracted from the state for each iteration.
data MonitorBatch a = MonitorBatch
  { mbName :: String,
    mbHandle :: Maybe Handle,
    mbParams :: [MonitorParameterBatch a],
    mbSize :: Int
  }

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

mbOpen :: String -> Bool -> MonitorBatch a -> IO (MonitorBatch a)
mbOpen n frc m = do
  let mfn = n <> mbName m <> ".batch"
  h <- open' mfn frc
  hSetBuffering h LineBuffering
  return $ m {mbHandle = Just h}

mbAppend :: String -> MonitorBatch a -> IO (MonitorBatch a)
mbAppend n m = do
  let fn = n <> mbName m <> ".batch"
  fe <- doesFileExist fn
  if fe
    then do
      h <- openFile fn AppendMode
      hSetBuffering h LineBuffering
      return $ m {mbHandle = Just h}
    else error $ "mbAppend: Monitor file does not exist: " ++ fn ++ "."

mbHeader :: MonitorBatch a -> IO ()
mbHeader m = case mbHandle m of
  Nothing ->
    error $
      "mbHeader: No handle available for batch monitor with name "
        <> mbName m
        <> "."
  Just h ->
    BL.hPutStrLn h $
      mfRenderRow $
        ["Iteration", "Mean log-Prior", "Mean log-Likelihood", "Mean log-Posterior"]
          ++ [BL.pack $ mbpName mbp | mbp <- mbParams m]

mean :: [Log Double] -> Log Double
mean xs = sum xs / fromIntegral (length xs)

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
      BL.hPutStrLn h $
        mfRenderRow $
          BL.pack (show i) :
          renderLog mlps :
          renderLog mlls :
          renderLog mlos :
            [BB.toLazyByteString $ mbpFunc mbp (map state t) | mbp <- mbParams m]
  where
    t = takeItems (mbSize m) t'
    lps = map prior t
    lls = map likelihood t
    los = zipWith (*) lps lls
    mlps = mean lps
    mlls = mean lls
    mlos = mean los

mbClose :: MonitorBatch a -> IO ()
mbClose m = case mbHandle m of
  Just h -> hClose h
  Nothing -> error $ "mfClose: File was not opened for batch monitor: " <> mbName m <> "."

-- | Open the files associated with the 'Monitor'.
mOpen :: String -> Bool -> Monitor a -> IO (Monitor a)
mOpen n frc (Monitor s fs bs) = do
  fs' <- mapM (mfOpen n frc) fs
  mapM_ mfHeader fs'
  bs' <- mapM (mbOpen n frc) bs
  mapM_ mbHeader bs'
  hSetBuffering stdout LineBuffering
  return $ Monitor s fs' bs'

-- | Open the files associated with the 'Monitor' in append mode.
mAppend :: String -> Monitor a -> IO (Monitor a)
mAppend n (Monitor s fs bs) = do
  fs' <- mapM (mfAppend n) fs
  bs' <- mapM (mbAppend n) bs
  return $ Monitor s fs' bs'

-- | Execute monitors; print status information to files and return text to be
-- printed to standard output and log file.
mExec ::
  -- | Verbosity
  Verbosity ->
  -- | Iteration.
  Int ->
  -- | Starting state.
  Int ->
  -- | Starting time.
  UTCTime ->
  -- | Trace of Markov chain.
  Trace a ->
  -- | Total number of iterations; to calculate ETA.
  Int ->
  -- | The monitor.
  Monitor a ->
  IO (Maybe BL.ByteString)
mExec v i ss st xs j (Monitor s fs bs) = do
  mapM_ (mfExec i $ headT xs) fs
  mapM_ (mbExec i xs) bs
  if v == Quiet
    then return Nothing
    else msExec i (headT xs) ss st j s

-- | Close the files associated with the 'Monitor'.
mClose :: Monitor a -> IO (Monitor a)
mClose m@(Monitor _ fms bms) = do
  mapM_ mfClose fms
  mapM_ mbClose bms
  let fms' = map (\fm -> fm {mfHandle = Nothing}) fms
  let bms' = map (\bm -> bm {mbHandle = Nothing}) bms
  return m {mFiles = fms', mBatches = bms'}
