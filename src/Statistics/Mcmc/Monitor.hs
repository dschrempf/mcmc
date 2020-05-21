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

module Statistics.Mcmc.Monitor
  ( Out
  , Monitor
  , monitorFile
  , monitorScreen
  , Control (..)
  , cOpen
  , cLog
  , cClose
  ) where

import qualified Data.ByteString.Char8 as B
import Data.ByteString (ByteString)
import System.IO

-- | Either log to file or handle.
data Out = OutFile FilePath | OutHandle Handle

outOpen :: Out -> IO Out
outOpen (OutFile   f) = OutHandle <$> openFile f WriteMode
outOpen (OutHandle h) = error $ "outOpen: Cannot open a handle: " <> show h <> "."

out :: Out -> ByteString -> IO ()
out (OutHandle h) = B.hPutStrLn h
out (OutFile   f) = error $ "out: No handle provided for file " <> f <> "."

-- TODO: Do not close stdout!
outClose :: Out -> IO ()
outClose (OutHandle h) = hClose h
outClose (OutFile   f) = error $ "outClose: Cannot close a file: " <> f <> "."

-- | Monitor a variable of the state space.
data Monitor a = Monitor
  {
    mOut  :: Out               -- ^ Log to file or handle ('stdout' can be used).
  , mShow :: a -> ByteString   -- ^ Instruction about what to log.
  , mFreq :: Int               -- ^ Logging frequency.
  }

mOpen :: Monitor a -> IO (Monitor a)
mOpen m = do
  h <- outOpen (mOut m)
  return $ m { mOut = h }

-- Execute monitor.
mLog :: Int -> a -> Monitor a -> IO ()
mLog i x m | i `mod` mFreq m /= 0 = return ()
           | otherwise = out (mOut m) (mShow m x)

mClose :: Monitor a -> IO ()
mClose m = do
  h <- outClose (mOut m)
  return ()

-- TODO: Write functions such as (open, print, close).

-- | Log to file.
monitorFile :: FilePath -> (a -> ByteString) -> Int -> Monitor a
monitorFile fp = Monitor (OutFile fp)

-- | Log variable to screen.
monitorScreen :: (a -> ByteString) -> Int -> Monitor a
monitorScreen = Monitor (OutHandle stdout)

-- TODO: Create monitors for specific types or type classes, such as Tree, or 'Num'.

-- | List of monitors.
newtype Control a = Control { fromControl :: [Monitor a] }

instance Semigroup (Control a) where
  (Control l) <> (Control r) = Control (l <> r)

instance Monoid (Control a)where
  mempty = Control []

-- | Open the files associated with the 'Control'.
cOpen :: Control a -> IO (Control a)
cOpen c = Control <$> mapM mOpen (fromControl c)

-- | Print logs for given iteration.
cLog :: Int -> a -> Control a -> IO ()
cLog i x = mapM_ (mLog i x) . fromControl

-- | Close the files associated with the 'Control'.
cClose :: Control a -> IO ()
cClose c = mapM_ mClose (fromControl c)
