{-# LANGUAGE RankNTypes #-}

-- TODO: Allow execution of moves in order of appearance in the cycle.

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
--
-- A short introduction to update mechanisms using the Metropolis-Hastings
-- algorithm (see Geyer, C. J., 2011; Introduction to Markov Chain Monte Carlo. In
-- Handbook of Markov Chain Monte Carlo (pp. 45), Chapman \& Hall/CRC).
--
-- A 'Move' is an instruction about how to advance a given Markov chain so that it
-- possibly reaches a new state. That is, 'Move's specify how the chain traverses
-- the state space. As far as this MCMC library is concerned, 'Move's are
-- /elementary updates/ in that they cannot be decomposed into smaller updates.
--
-- 'Move's can be combined to form composite updates, a technique often referred to
-- as /composition/. On the other hand, /mixing/ (used in the sense of mixture
-- models) is the random choice of a 'Move' (or a composition of 'Move's) from a
-- given set.
--
-- The __composition__ and __mixture__ of 'Move's allows specification of nearly
-- all MCMC algorithms involving a single chain (i.e., population methods such as
-- particle filters are excluded). In particular, Gibbs samplers of all sorts can
-- be specified using this procedure.
--
-- This library enables composition and mixture of 'Move's via the 'Cycle' data
-- type. Essentially, a 'Cycle' is a set of 'Move's. The chain advances after the
-- completion of each 'Cycle', which is called an __iteration__, and the iteration
-- counter is increased by one.
--
-- The 'Move's in a 'Cycle' can be executed in the given order or in a random
-- sequence which allows, for example, specification of a fixed scan Gibbs sampler,
-- or a random sequence scan Gibbs sampler, respectively.
--
-- Note that it is of utter importance that the given 'Cycle' enables traversal of
-- the complete state space. Otherwise, the Markov chain will not converge to the
-- correct stationary posterior distribution.
module Mcmc.Move
  ( -- * Move
    Move (..),
    MoveSimple (..),
    Tuner (tParam, tFunc),
    tuner,
    tune,
    autotune,

    -- * Cycle
    Cycle (fromCycle),
    fromList,
    autotuneC,
    summarizeCycle,

    -- * Acceptance
    Acceptance (..),
    emptyA,
    prependA,
    resetA,
    acceptanceRatios,
  )
where

import Data.Aeson
import Data.Function
import Data.List
import qualified Data.Map.Strict as M
import Data.Map.Strict (Map)
import Data.Maybe
import Numeric.Log hiding (sum)
import System.Random.MWC
import Text.Printf

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
  show = mvName

instance Eq (Move a) where
  m == n = mvName m == mvName n

instance Ord (Move a) where
  compare = compare `on` mvName

-- XXX: One could also use a different type for 'mvSample', so that
-- 'mvLogDensity' can be avoided. In detail,
--
-- @
--   mvSample :: a -> GenIO -> IO (a, Log Double, Log, Double)
-- @
--
-- where the log densities describe the probability of going there and back.
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
    mvSample ::
      a ->
      GenIO ->
      IO a,
    -- | The log-density of going from one state to another.
    mvLogDensity ::
      a ->
      a ->
      Log Double
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
tune dt tm
  | dt <= 0 = error $ "tune: Tuning parameter not positive: " <> show dt <> "."
  | otherwise = do
    (Tuner t f) <- mvTune tm
    let t' = t * dt
    return $ tm {mvSimple = f t', mvTune = Just $ Tuner t' f}

ratioOpt :: Double
ratioOpt = 0.44

-- | For a given acceptance ratio, auto tune the 'Move'. For now, a 'Move' is
-- enlarged when the acceptance ratio is above 0.44, and shrunk otherwise.
-- Return 'Nothing' if 'Move' is not tuneable.
--
-- XXX: The desired acceptance ratio 0.44 is optimal for one-dimensional
-- 'Move's; one could also store the affected number of dimensions with the
-- 'Move' and tune towards an acceptance ratio accounting for the number of
-- dimensions.
autotune :: Double -> Move a -> Maybe (Move a)
autotune a = tune (a / ratioOpt)

-- | In brief, a 'Cycle' is a list of moves. The state of the Markov chain will
-- be logged only after each 'Cycle', and the iteration counter will be
-- increased by one. __Moves must have unique names__, so that they can be
-- identified.
--
-- 'Move's are replicated according to their weights and executed in random
-- order. That is, they are not executed in the order they appear in the
-- 'Cycle'. However, if a 'Move' has weight @w@, it is executed exactly @w@
-- times per iteration.
newtype Cycle a = Cycle {fromCycle :: [Move a]}

-- | Create a 'Cycle' from a list of 'Move's.
fromList :: [Move a] -> Cycle a
fromList [] =
  error "fromList: Received an empty list but cannot create an empty Cycle."
fromList xs =
  if length (nub nms) == length nms
    then Cycle xs
    else error "fromList: Moves don't have unique names."
  where
    nms = map mvName xs

-- | Tune the 'Move's in the 'Cycle'. Tuning has no effect on 'Move's that
-- cannot be tuned. See 'autotune'.
autotuneC :: Int -> Acceptance (Move a) -> Cycle a -> Cycle a
autotuneC n a = Cycle . map tuneF . fromCycle
  where
    ars = acceptanceRatios n a
    tuneF m = fromMaybe m (autotune (ars M.! m) m)

left :: Int -> String -> String
left n s = take n s ++ replicate (n - l) ' ' where l = length s

right :: Int -> String -> String
right n s = replicate (n - l) ' ' ++ take n s where l = length s

renderRow :: String -> String -> String -> String -> String
renderRow name weight acceptRatio tuneParam = nB ++ wB ++ rB ++ tB
  where
    nB = left 30 name
    wB = right 8 weight
    rB = right 20 acceptRatio
    tB = right 20 tuneParam

moveHeader :: String
moveHeader =
  renderRow "Move name" "Weight" "Acceptance ratio" "Tuning parameter"

summarizeMove :: Move a -> Maybe Double -> String
summarizeMove m r = renderRow name weight acceptRatio tuneParamStr
  where
    name = mvName m
    weight = show $ mvWeight m
    acceptRatio = maybe "" (printf "%.3f") r
    tuneParamStr = maybe "" (printf "%.3f") (tParam <$> mvTune m)

-- TODO: Move this to Mcmc?

-- | Summarize the 'Move's in the 'Cycle'. Also report acceptance ratios for the
-- given number of last iterations.
summarizeCycle :: Maybe (Int, Acceptance (Move a)) -> Cycle a -> String
summarizeCycle acc c =
  unlines $
    [ "Summary of move(s) in cycle:",
      replicate (length moveHeader) '─',
      moveHeader,
      replicate (length moveHeader) '─'
    ]
      ++ [summarizeMove m (ar m) | m <- mvs]
      ++ [ replicate (length moveHeader) '─',
           show mpi ++ " move(s) per iteration." ++ arStr
         ]
  where
    mvs = fromCycle c
    mpi = sum $ map mvWeight mvs
    arStr = case acc of
      Nothing -> ""
      Just (n, _) ->
        " Acceptance ratio(s) calculated over " ++ show n ++ " iterations."
    ars = case acc of
      Nothing -> M.empty
      Just (n, a) -> acceptanceRatios n a
    ar m = ars M.!? m

-- XXX: I am not too happy about this data type but cannot think of a better
-- solution.

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
prependA :: (Ord k, Show k) => k -> Bool -> Acceptance k -> Acceptance k
-- Unsafe; faster.
prependA k v (Acceptance m) = Acceptance $ M.adjust (v :) k m
-- -- Safe; slower.
-- prependA k v (Acceptance m) | k `M.member` m = Acceptance $ M.adjust (v:) k m
--                             | otherwise = error msg
--   where msg = "prependA: Can not add acceptance value for key: " <> show k <> "."
{-# INLINEABLE prependA #-}

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
