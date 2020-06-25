{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE RankNTypes #-}

-- TODO: Moves and monitors both use lenses and names for what they move and
-- monitor. Should a data structure be used combining the lens and the name, so
-- that things are cohesive?

-- TODO: Moves on simplices: SimplexElementScale (?).

-- TODO: Moves on tree branch lengths.
-- - Slide a node on the tree.
-- - Scale a tree.

-- TODO: Moves on tree topologies.
-- - NNI
-- - Narrow (what is this, see RevBayes)
-- - FNPR (dito)

-- TODO: Bactrian moves; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/.
--
-- slideBactrian
--
-- scaleBactrian

-- TODO: Try to use abstract data types (Acceptance, Move, MoveSimple).

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
    autotune,

    -- * Cycle
    MoveOrder (..),
    Cycle (ccMoves),
    fromList,
    setOrder,
    getNCycles,
    tuneCycle,
    autotuneCycle,
    summarizeCycle,

    -- * Acceptance
    Acceptance (..),
    emptyA,
    pushA,
    resetA,
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
    mvTune :: Maybe (Tuner a)
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

-- | Tune the acceptance ratio of a 'Move'; see 'tune', or 'autotune'.
data Tuner a = Tuner
  { tParam :: Double,
    tFunc :: Double -> MoveSimple a
  }

-- | Create a 'Tuner'. The tuning function accepts a tuning parameter, and
-- returns a corresponding 'MoveSimple'. The larger the tuning parameter, the
-- larger the 'Move', and vice versa.
tuner :: (Double -> MoveSimple a) -> Tuner a
tuner = Tuner 1.0

-- | Tune a 'Move'. Return 'Nothing' if 'Move' is not tuneable. If the parameter
--   @dt@ is larger than 1.0, the 'Move' is enlarged, if @0<dt<1.0@, it is
--   shrunk. Negative tuning parameters are not allowed.
tune :: Double -> Move a -> Maybe (Move a)
tune dt m
  | dt <= 0 = error $ "tune: Tuning parameter not positive: " <> show dt <> "."
  | otherwise = do
    (Tuner t f) <- mvTune m
    let t' = t * dt
    return $ m {mvSimple = f t', mvTune = Just $ Tuner t' f}

-- XXX: The desired acceptance ratio 0.44 is optimal for one-dimensional
-- 'Move's; one could also store the affected number of dimensions with the
-- 'Move' and tune towards an acceptance ratio accounting for the number of
-- dimensions.
ratioOpt :: Double
ratioOpt = 0.44

-- | For a given acceptance ratio, auto tune the 'Move'. For now, a 'Move' is
-- enlarged when the acceptance ratio is above 0.44, and shrunk otherwise.
-- Return 'Nothing' if 'Move' is not tuneable.
autotune :: Double -> Move a -> Maybe (Move a)
autotune a = tune (a / ratioOpt)

-- | Define the order in which 'Move's are executed in a 'Cycle'. Technically,
-- this is not an order, since the total number of 'Move's per 'Cycle' may be
-- different (e.g., compare 'RandomO' and 'RandomReversibleO').
data MoveOrder
  = -- | Shuffle the 'Move's in the 'Cycle'. The 'Move's are replicated
    -- according to their weights and executed in random order. If a 'Move' has
    -- weight @w@, it is executed exactly @w@ times per iteration.
    RandomO
  | -- | The 'Move's are executed sequentially, in the order they appear in the
    -- 'Cycle'. 'Move's with weight @w>1@ are repeated immediately @w@ times
    -- (and not appended to the end of the list).
    SequentialO
  | -- | Similar to 'RandomOrder'. However, a reversed copy of the list of
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

instance Default MoveOrder where def = RandomO

-- | In brief, a 'Cycle' is a list of moves. The state of the Markov chain will
-- be logged only after all 'Move's in the 'Cycle' have been completed, and the
-- iteration counter will be increased by one. The order in which the 'Move's
-- are executed is specified by 'MoveOrder'. The default is 'RandomOrder'.
--
-- __Moves must have unique names__, so that they can be identified.
data Cycle a = Cycle
  { ccMoves :: [Move a],
    ccOrder :: MoveOrder
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
setOrder :: MoveOrder -> Cycle a -> Cycle a
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

-- | Tune some 'Move's in the 'Cycle'. See 'tune'. Moves in the map, but not in
-- the cycle are ignored!
tuneCycle :: Map (Move a) Double -> Cycle a -> Cycle a
tuneCycle m c = c {ccMoves = map tuneF $ ccMoves c}
  where
    tuneF mv = case m M.!? mv of
      Nothing -> mv
      Just x -> fromMaybe mv (tune x mv)

-- | Caculate acceptance ratios for the given number of last iterations. Auto
-- tune the 'Move's in the 'Cycle' with the calculated acceptance ratios. See
-- 'autotune'.
autotuneCycle :: Int -> Acceptance (Move a) -> Cycle a -> Cycle a
autotuneCycle n a = tuneCycle (M.map (/ratioOpt) $ acceptanceRatios n a)

renderRow :: Text -> Text -> Text -> Text -> Text
renderRow name weight acceptRatio tuneParam = "   " <> nB <> wB <> rB <> tB
  where
    nB = T.justifyLeft 30 ' ' name
    wB = T.justifyRight 8 ' ' weight
    rB = T.justifyRight 20 ' ' acceptRatio
    tB = T.justifyRight 20 ' ' tuneParam

moveHeader :: Text
moveHeader =
  renderRow "Move name" "Weight" "Acceptance ratio" "Tuning parameter"

summarizeMove :: Move a -> Maybe Double -> Text
summarizeMove m r = renderRow (T.pack name) weight acceptRatio tuneParamStr
  where
    name = mvName m
    weight = B.toLazyText $ B.decimal $ mvWeight m
    acceptRatio = B.toLazyText $ maybe "" (B.formatRealFloat B.Fixed (Just 3)) r
    tuneParamStr = B.toLazyText $ maybe "" (B.formatRealFloat B.Fixed (Just 3)) (tParam <$> mvTune m)

-- | Summarize the 'Move's in the 'Cycle'. Also report acceptance ratios for the
-- given number of last iterations.
summarizeCycle :: Maybe (Int, Acceptance (Move a)) -> Cycle a -> Text
summarizeCycle acc c =
  T.unlines $
    [ "-- Summary of move(s) in cycle.",
      -- T.replicate (T.length moveHeader) "─",
      moveHeader,
      "   " <> T.replicate (T.length moveHeader - 3) "─"
    ]
      ++ [summarizeMove m (ar m) | m <- mvs]
      ++ [ "   " <> T.replicate (T.length moveHeader - 3) "─",
           B.toLazyText $
             B.fromLazyText "-- "
               <> B.decimal mpi
               <> B.fromString " move(s) per iteration."
               <> arStr
         ]
  where
    mvs = ccMoves c
    mpi = sum $ map mvWeight mvs
    arStr = case acc of
      Nothing -> ""
      Just (n, _) ->
        " Acceptance ratio(s) calculated over " <> B.decimal n <> " iterations."
    ars = case acc of
      Nothing -> M.empty
      Just (n, a) -> acceptanceRatios n a
    ar m = ars M.!? m

-- | For each key @k@, store the list of accepted (True) and rejected (False)
-- proposals. For reasons of efficiency, the lists are stored in reverse order;
-- latest first.
newtype Acceptance k = Acceptance {fromAcceptance :: Map k [Bool]}

instance ToJSONKey k => ToJSON (Acceptance k) where
  toJSON (Acceptance m) = toJSON m
  toEncoding (Acceptance m) = toEncoding m

instance (Ord k, FromJSONKey k) => FromJSON (Acceptance k) where
  parseJSON v = Acceptance <$> parseJSON v

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
emptyA :: Ord k => [k] -> Acceptance k
emptyA ks = Acceptance $ M.fromList [(k, []) | k <- ks]

-- | For key @k@, prepend an accepted (True) or rejected (False) proposal.
pushA :: (Ord k, Show k) => k -> Bool -> Acceptance k -> Acceptance k
-- Unsafe; faster.
pushA k v (Acceptance m) = Acceptance $ M.adjust (v :) k m
-- -- Safe; slower.
-- prependA k v (Acceptance m) | k `M.member` m = Acceptance $ M.adjust (v:) k m
--                             | otherwise = error msg
--   where msg = "prependA: Can not add acceptance value for key: " <> show k <> "."
{-# INLINEABLE pushA #-}

-- | Reset acceptance storage.
resetA :: Acceptance k -> Acceptance k
resetA = Acceptance . M.map (const []) . fromAcceptance

ratio :: Int -> [Bool] -> Double
ratio n xs = fromIntegral (length ts) / fromIntegral n
  where
    ts = filter (== True) xs

-- | Acceptance ratios averaged over the given number of last iterations. If
-- less than @n@ iterations are available, only those are used.
acceptanceRatios :: Int -> Acceptance k -> Map k Double
acceptanceRatios n (Acceptance m) = M.map (ratio n . take n) m
