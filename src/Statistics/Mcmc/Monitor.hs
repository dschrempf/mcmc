{-# LANGUAGE OverloadedStrings #-}

{- |
Module      :  Statistics.Mcmc.Monitor
Description :  Monitor a Markov chain
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 21 14:35:11 2020.

-}

-- TODO: Improvements necessary.
-- Something like: one monitor prints a set of parameters.

-- TODO: Use builders.

module Statistics.Mcmc.Monitor
  ( Out
  , Monitor
  , monitorFile
  , monitorStdOut
  , Monitors (..)
  , msOpen
  , msExec
  , msClose
  ) where

import qualified Data.ByteString.Char8 as B
import Data.ByteString (ByteString)
import System.IO

-- | Either log to file, or to a handle. Before the analysis, the given file
-- will be opened using 'WriteMode' to yield a handle; the handle will be closed
-- at the end. If a handle is given, this handle will be used for output. In
-- this case, the handle will not be closed at the end.
data Out = OFile FilePath | OHandle Handle
  deriving (Show)

getHandle :: Out -> IO Handle
getHandle (OFile   f) = openFile f WriteMode
getHandle (OHandle h) = pure h

-- | Monitor a variable of the state space.
data Monitor a = Monitor
  {
    mOut    :: Out               -- ^ Log to file or standard output.
  , mHandle :: Maybe Handle      -- ^ Handle to use.
  , mShow   :: a -> ByteString   -- ^ Instruction about what to log.
  , mFreq   :: Int               -- ^ Logging period.
  }

mOpen :: Monitor a -> IO (Monitor a)
mOpen m = do
  h <- getHandle (mOut m)
  return $ m { mHandle = Just h }

mExec :: Int -> a -> Monitor a -> IO ()
mExec i x m | i `mod` mFreq m /= 0 = return ()
            | otherwise = case mHandle m of
                Just h  -> B.hPutStrLn h $ B.pack (show i) <> "\t" <> mShow m x
                Nothing -> error $ "mLog: No handle available for monitor " <> show (mOut m) <> "."

mClose :: Monitor a -> IO ()
mClose m = case (mOut m, mHandle m) of
  (OFile   _, Just h ) -> hClose h
  (OFile   f, Nothing) -> error $ "mClose: File was not opened: " <> f <> "."
  (OHandle _, _      ) -> pure ()

-- | Log to file.
monitorFile
  :: FilePath          -- ^ File path; file will be overwritten!
  -> (a -> ByteString) -- ^ Instruction about what to log.
  -> Int               -- ^ Logging period.
  -> Monitor a
monitorFile fp = Monitor (OFile fp) Nothing

-- | Log to 'stdout'; see 'monitorFile'.
monitorStdOut :: (a -> ByteString) -> Int -> Monitor a
monitorStdOut = Monitor (OHandle stdout) Nothing

-- TODO: Create monitors for specific types or type classes, such as Tree, or 'Num'.

-- | List of monitors.
newtype Monitors a = Monitors { fromMonitors :: [Monitor a] }

instance Semigroup (Monitors a) where
  (Monitors l) <> (Monitors r) = Monitors (l <> r)

instance Monoid (Monitors a) where
  mempty = Monitors []

-- | Open the files associated with the 'Monitors'.
msOpen :: Monitors a -> IO (Monitors a)
msOpen c = Monitors <$> mapM mOpen (fromMonitors c)

-- | Print logs for given iteration.
msExec :: Int -> a -> Monitors a -> IO ()
msExec i x = mapM_ (mExec i x) . fromMonitors

-- | Close the files associated with the 'Monitors'.
msClose :: Monitors a -> IO ()
msClose c = mapM_ mClose (fromMonitors c)
