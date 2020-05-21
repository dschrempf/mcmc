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

-- | Either log to file (which will be opened using 'WriteMode' to yield a
-- handle), or standard output.
data Out = File FilePath | StdOut
  deriving (Show)

getHandle :: Out -> IO Handle
getHandle (File f) = openFile f WriteMode
getHandle StdOut   = pure stdout

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

mLog :: Int -> a -> Monitor a -> IO ()
mLog i x m | i `mod` mFreq m /= 0 = return ()
           | otherwise = case mHandle m of
               Just h  -> B.hPutStrLn h (mShow m x)
               Nothing -> error $ "mLog: No handle available for monitor " <> show (mOut m) <> "."

mClose :: Monitor a -> IO ()
mClose m = case (mOut m, mHandle m) of
  (File _, Just h ) -> hClose h
  (StdOut, _      ) -> pure ()
  (File f, Nothing) -> error $ "mClose: File was not opened: " <> f <> "."

-- | Log to file.
monitorFile
  :: FilePath          -- ^ File path; file will be overwritten!
  -> (a -> ByteString) -- ^ Instruction about what to log.
  -> Int               -- ^ Logging period.
  -> Monitor a
monitorFile fp = Monitor (File fp) Nothing

-- | Log to screen; see 'monitorFile'.
monitorScreen :: (a -> ByteString) -> Int -> Monitor a
monitorScreen = Monitor StdOut Nothing

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
