{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- TODO: Moves on simplices: SimplexElementScale (?).

-- TODO: Moves on tree branch lengths.
-- - Slide a node on the tree.
-- - Scale a tree.

-- TODO: Moves on tree topologies.
-- - NNI
-- - Narrow (what is this, see RevBayes)
-- - FNPR (dito)

-- |
-- Module      :  Mcmc.Move
-- Description :  Moves and cycles
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 13:42:53 2020.
module Mcmc.Move
  ( -- * Move
    Move (..),
    MoveSimple (..),
    Tuner (tParam, tFunc),
    tuner,
    tune,

    -- * Cycle
    Order (..),
    Cycle (ccMoves),
    fromList,
    setOrder,
    getNCycles,
    tuneCycle,
    autotuneCycle,
    summarizeCycle,

    -- * Acceptance
    Acceptance (fromAcceptance),
    emptyA,
    pushA,
    resetA,
    transformKeysA,
    acceptanceRatios,
  )
where

import Data.Aeson
import Data.Default
import Data.Function
import Data.List
import qualified Data.Map.Strict as M
import Data.Map.Strict (Map)
import Data.Maybe
import qualified Data.Text.Lazy as T
import Data.Text.Lazy (Text)
import qualified Data.Text.Lazy.Builder as B
import qualified Data.Text.Lazy.Builder.Int as B
import qualified Data.Text.Lazy.Builder.RealFloat as B
import Lens.Micro
import Mcmc.Tools.Shuffle
import Numeric.Log hiding (sum)
import System.Random.MWC

-- | A 'Move' is an instruction about how the Markov chain will traverse the
-- state space @a@. Essentially, it is a probability density conditioned on the
-- current state.
--
-- A 'Move' may be tuneable in that it contains information about how to enlarge
-- or shrink the step size to tune the acceptance ratio.
data Move a = Move
  { -- | Name (no moves with the same name are allowed in a 'Cycle').
    mvName :: String,
    -- | The weight determines how often a 'Move' is executed per iteration of
    -- the Markov chain.
    mvWeight :: Int,
    -- | Simple move without tuning information.
    mvSimple :: MoveSimple a,
    -- | Tuning is disabled if set to 'Nothing'.
    mvTuner :: Maybe (Tuner a)
  }

instance Show (Move a) where
  show m = show $ mvName m

instance Eq (Move a) where
  m == n = mvName m == mvName n

instance Ord (Move a) where
  compare = compare `on` mvName

-- One could also use a different type for 'mvSample', so that 'mvDensity' can
-- be avoided. In detail,
--
-- @
--   mvSample :: a -> GenIO -> IO (a, Log Double, Log, Double)
-- @
--
-- where the densities describe the probability of going there and back.
-- However, we may need more information about the move for other MCMC samplers
-- different from Metropolis-Hastings.

-- | Simple move without tuning information.
--
-- In order to calculate the Metropolis-Hastings ratio, we need to know the
-- probability (density) of jumping forth, and the probability (density) of
-- jumping back.
data MoveSimple a = MoveSimple
  { -- | Instruction about randomly moving from the current state to a new
    -- state, given some source of randomness.
    mvSample :: a -> GenIO -> IO a,
    -- | The density of going from one state to another. Set to 'Nothing' for
    -- symmetric moves.
    mvDensity :: Maybe (a -> a -> Log Double)
  }

-- | Tune the acceptance ratio of a 'Move'; see 'tune', or 'autotuneCycle'.
data Tuner a = Tuner
  { tParam :: Double,
    tFunc :: Double -> MoveSimple a
  }

-- | Create a 'Tuner'. The tuning function accepts a tuning parameter, and
-- returns a corresponding 'MoveSimple'. The larger the tuning parameter, the
-- larger the 'Move', and vice versa.
tuner :: (Double -> MoveSimple a) -> Tuner a
tuner = Tuner 1.0

-- Minimal tuning parameter; subject to change.
tuningParamMin :: Double
tuningParamMin = 1e-12

-- | Tune a 'Move'. Return 'Nothing' if 'Move' is not tuneable. If the parameter
--   @dt@ is larger than 1.0, the 'Move' is enlarged, if @0<dt<1.0@, it is
--   shrunk. Negative tuning parameters are not allowed.
tune :: Double -> Move a -> Maybe (Move a)
tune dt m
  | dt <= 0 = error $ "tune: Tuning parameter not positive: " <> show dt <> "."
  | otherwise = do
    (Tuner t f) <- mvTuner m
    -- Ensure that the tuning parameter is not too small.
    let t' = max tuningParamMin (t * dt)
    return $ m {mvSimple = f t', mvTuner = Just $ Tuner t' f}

-- XXX: The desired acceptance ratio 0.44 is optimal for one-dimensional
-- 'Move's; one could also store the affected number of dimensions with the
-- 'Move' and tune towards an acceptance ratio accounting for the number of
-- dimensions.
ratioOpt :: Double
ratioOpt = 0.44

-- | Define the order in which 'Move's are executed in a 'Cycle'. The total
-- number of 'Move's per 'Cycle' may differ between 'Order's (e.g., compare
-- 'RandomO' and 'RandomReversibleO').
data Order
  = -- | Shuffle the 'Move's in the 'Cycle'. The 'Move's are replicated
    -- according to their weights and executed in random order. If a 'Move' has
    -- weight @w@, it is executed exactly @w@ times per iteration.
    RandomO
  | -- | The 'Move's are executed sequentially, in the order they appear in the
    -- 'Cycle'. 'Move's with weight @w>1@ are repeated immediately @w@ times
    -- (and not appended to the end of the list).
    SequentialO
  | -- | Similar to 'RandomO'. However, a reversed copy of the list of
    --  shuffled 'Move's is appended such that the resulting Markov chain is
    --  reversible.
    --  Note: the total number of 'Move's executed per cycle is twice the number
    --  of 'RandomO'.
    RandomReversibleO
  | -- | Similar to 'SequentialO'. However, a reversed copy of the list of
    -- sequentially ordered 'Move's is appended such that the resulting Markov
    -- chain is reversible.
    SequentialReversibleO
  deriving (Eq, Show)

instance Default Order where def = RandomO

-- | In brief, a 'Cycle' is a list of moves. The state of the Markov chain will
-- be logged only after all 'Move's in the 'Cycle' have been completed, and the
-- iteration counter will be increased by one. The order in which the 'Move's
-- are executed is specified by 'Order'. The default is 'RandomO'.
--
-- __Moves must have unique names__, so that they can be identified.
data Cycle a = Cycle
  { ccMoves :: [Move a],
    ccOrder :: Order
  }

-- | Create a 'Cycle' from a list of 'Move's.
fromList :: [Move a] -> Cycle a
fromList [] =
  error "fromList: Received an empty list but cannot create an empty Cycle."
fromList xs =
  if length (nub nms) == length nms
    then Cycle xs def
    else error "fromList: Moves don't have unique names."
  where
    nms = map mvName xs

-- | Set the order of 'Move's in a 'Cycle'.
setOrder :: Order -> Cycle a -> Cycle a
setOrder o c = c {ccOrder = o}

-- | Replicate 'Move's according to their weights and possibly shuffle them.
getNCycles :: Cycle a -> Int -> GenIO -> IO [[Move a]]
getNCycles (Cycle xs o) n g = case o of
  RandomO -> shuffleN mvs n g
  SequentialO -> return $ replicate n mvs
  RandomReversibleO -> do
    mvsRs <- shuffleN mvs n g
    return [mvsR ++ reverse mvsR | mvsR <- mvsRs]
  SequentialReversibleO -> return $ replicate n $ mvs ++ reverse mvs
  where
    !mvs = concat [replicate (mvWeight m) m | m <- xs]

-- | Tune 'Move's in the 'Cycle'. See 'tune'.
tuneCycle :: Map (Move a) Double -> Cycle a -> Cycle a
tuneCycle m c =
  if sort (M.keys m) == sort mvs
    then c {ccMoves = map tuneF mvs}
    else error "tuneCycle: Map contains moves that are not in the cycle."
  where
    mvs = ccMoves c
    tuneF mv = case m M.!? mv of
      Nothing -> mv
      Just x -> fromMaybe mv (tune x mv)

-- | Calculate acceptance ratios and auto tune the 'Move's in the 'Cycle'. For
-- now, a 'Move' is enlarged when the acceptance ratio is above 0.44, and
-- shrunk otherwise. Do not change 'Move's that are not tuneable.
autotuneCycle :: Acceptance (Move a) -> Cycle a -> Cycle a
autotuneCycle a = tuneCycle (M.map (\x -> exp $ x - ratioOpt) $ acceptanceRatios a)

renderRow :: Text -> Text -> Text -> Text -> Text -> Text -> Text
renderRow name weight nAccept nReject acceptRatio tuneParam = "   " <> nm <> wt <> na <> nr <> ra <> tp
  where
    nm = T.justifyLeft 30 ' ' name
    wt = T.justifyRight 8 ' ' weight
    na = T.justifyRight 15 ' ' nAccept
    nr = T.justifyRight 15 ' ' nReject
    ra = T.justifyRight 15 ' ' acceptRatio
    tp = T.justifyRight 20 ' ' tuneParam

moveHeader :: Text
moveHeader =
  renderRow "Move name" "Weight" "Accepted" "Rejected" "Ratio" "Tuning parameter"

summarizeMove :: Move a -> Maybe (Int, Int, Double) -> Text
summarizeMove m r = renderRow (T.pack name) weight nAccept nReject acceptRatio tuneParamStr
  where
    name = mvName m
    weight = B.toLazyText $ B.decimal $ mvWeight m
    nAccept = B.toLazyText $ maybe "" (B.decimal . (^. _1)) r
    nReject = B.toLazyText $ maybe "" (B.decimal . (^. _2)) r
    acceptRatio = B.toLazyText $ maybe "" (B.formatRealFloat B.Fixed (Just 3) . (^. _3)) r
    tuneParamStr = B.toLazyText $ maybe "" (B.formatRealFloat B.Fixed (Just 3)) (tParam <$> mvTuner m)

-- | Summarize the 'Move's in the 'Cycle'. Also report acceptance ratios.
summarizeCycle :: Acceptance (Move a) -> Cycle a -> Text
summarizeCycle a c =
  T.unlines $
    [ "-- Summary of move(s) in cycle. " <> mpi <> " move(s) per iteration.",
      moveHeader,
      "   " <> T.replicate (T.length moveHeader - 3) "─"
    ]
      ++ [summarizeMove m (ar m) | m <- mvs]
      ++ [ "   " <> T.replicate (T.length moveHeader - 3) "─"]
  where
    mvs = ccMoves c
    mpi = B.toLazyText $ B.decimal $ sum $ map mvWeight mvs
    ar m = acceptanceRatio m a

-- | For each key @k@, store the number of accepted and rejected proposals.
newtype Acceptance k = Acceptance {fromAcceptance :: Map k (Int, Int)}

instance ToJSONKey k => ToJSON (Acceptance k) where
  toJSON (Acceptance m) = toJSON m
  toEncoding (Acceptance m) = toEncoding m

instance (Ord k, FromJSONKey k) => FromJSON (Acceptance k) where
  parseJSON v = Acceptance <$> parseJSON v

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
emptyA :: Ord k => [k] -> Acceptance k
emptyA ks = Acceptance $ M.fromList [(k, (0, 0)) | k <- ks]

-- | For key @k@, prepend an accepted (True) or rejected (False) proposal.
pushA :: (Ord k, Show k) => k -> Bool -> Acceptance k -> Acceptance k
pushA k True = Acceptance . M.adjust (\(a, r) -> (succ a, r)) k . fromAcceptance
pushA k False = Acceptance . M.adjust (\(a, r) -> (a, succ r)) k . fromAcceptance
{-# INLINEABLE pushA #-}

-- | Reset acceptance storage.
resetA :: Ord k => Acceptance k -> Acceptance k
resetA = emptyA . M.keys . fromAcceptance

transformKeys :: (Ord k1, Ord k2) => [k1] -> [k2] -> Map k1 v -> Map k2 v
transformKeys ks1 ks2 m = foldl' insrt M.empty $ zip ks1 ks2
  where
    insrt m' (k1, k2) = M.insert k2 (m M.! k1) m'

-- | Transform keys using the given lists. Keys not provided will not be present
-- in the new 'Acceptance' variable.
transformKeysA :: (Ord k1, Ord k2) => [k1] -> [k2] -> Acceptance k1 -> Acceptance k2
transformKeysA ks1 ks2 = Acceptance . transformKeys ks1 ks2 . fromAcceptance

-- | Acceptance counts and ratio for a specific move.
acceptanceRatio :: (Show k, Ord k) => k -> Acceptance k -> Maybe (Int, Int, Double)
acceptanceRatio k a = case fromAcceptance a M.!? k of
  Just (0, 0) -> Nothing
  Just (as, rs) -> Just (as, rs, fromIntegral as / fromIntegral (as + rs))
  Nothing -> error $ "acceptanceRatio: Key not found in map: " ++ show k ++ "."

-- | Acceptance ratios for all moves.
acceptanceRatios :: Acceptance k -> Map k Double
acceptanceRatios = M.map (\(as, rs) -> fromIntegral as / fromIntegral (as + rs)) . fromAcceptance
