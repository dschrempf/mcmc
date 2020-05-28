{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Statistics.Mcmc.Move
Description :  Moves and cycles
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed May 20 13:42:53 2020.

'Move's specify how the chain traverses the state space. A set of 'Move's is
called a 'Cycle'. The chain advances after the completion of each 'Cycle', and
the iteration counter is increased by one. It is important that the given
'Cycle' enables traversal of the complete state space. Otherwise, the Markov
chain will not converge to the correct stationary posterior distribution.

-}

-- TODO: Iteration -> cycle.

module Statistics.Mcmc.Move
  (
    -- * Move
    Move (..)
  , MoveSimple (..)
  , Tuner (tParam, tFunc)
  , tuner
  , tune
  , autotune
    -- * Cycle
  , Cycle (fromCycle)
  , addMove
  , fromList
  , mapCycle
  , summarizeCycle
  , summarizeCycleA
    -- * Acceptance
  , Acceptance
  , empty
  , prependA
  , resetA
  , acceptanceRatios
  ) where

import Data.Function
import qualified Data.Map.Strict as M
import Data.Map.Strict (Map)
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
  {
    -- | Name (no moves with the same name are allowed in a 'Cycle').
    mvName   :: String
    -- | The weight determines how often a 'Move' is executed per 'Cycle'.
  , mvWeight :: Int
    -- | Simple move without tuning information.
  , mvSimple :: MoveSimple a
    -- | Tuning is disabled if set to 'Nothing'.
  , mvTune   :: Maybe (Tuner a)
  }

instance Show (Move a) where
  show = mvName

instance Eq (Move a) where
  m == n =
    mvName m == mvName n

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
  {
    -- | Instruction about randomly moving from the current state to a new
    -- state, given some source of randomness.
    mvSample     :: a
                 -> GenIO
                 -> IO a
    -- | The log-density of going from one state to another.
  , mvLogDensity :: a
                 -> a
                 -> Log Double
  }

-- | Tune the acceptance ratio of a 'Move'; see 'tune', or 'autotune'.
data Tuner a = Tuner
  {
    tParam :: Double
  , tFunc  :: Double -> MoveSimple a
  }

-- | Create a 'Tuner'. The tuning function accepts a tuning parameter, and
-- returns a corresponding 'MoveSimple'. The larger the tuning parameter, the
-- larger the 'Move', and vice versa.
tuner :: (Double -> MoveSimple a) -> Tuner a
tuner = Tuner 1.0

-- | Tune a 'Move'. Return 'Nothing' if 'Move' is not tuneable. If the parameter
--   @dt@ is larger than 1.0, the 'Move' is enlarged, if @0<dt<1.0@, it is
--   shrunk. Negative tuning parameters are not allowed.
tune
  :: Double
  -> Move a
  -> Maybe (Move a)
tune dt tm
  | dt <= 0   = error $ "tune: Tuning parameter not positive: " <> show dt <> "."
  | otherwise = do
      (Tuner t f) <- mvTune tm
      let t' = t*dt
      return $ tm { mvSimple = f t'
                  , mvTune = Just $ Tuner t' f }

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
-- XXX: The tuning parameter is experimental; other distributions could be
-- tried.
autotune a = tune (a / ratioOpt)

-- | In brief, a 'Cycle' is a list of moves. The state of the Markov chain will
-- be logged only after each 'Cycle'. __Moves must have unique names__, so that
-- they can be identified.
--
-- 'Move's are replicated according to their weights and executed in random
-- order. That is, they are not executed in the order they appear in the
-- 'Cycle'.
newtype Cycle a = Cycle { fromCycle :: [Move a] }
-- XXX: A classical list is used, and not 'Data.List.NonEmpty' (keep things
-- simple).

-- | Add a 'Move' to the 'Cycle'. The name of the added 'Move' must be unique.
addMove :: Move a -> Cycle a -> Cycle a
addMove m c | m `notElem` fromCycle c = Cycle $ m : fromCycle c
            | otherwise = error msg
  where msg = "addMove: Move " <> mvName m <> " already exists in cycle."

-- | Create a 'Cycle' from a list of 'Move's.
fromList :: [Move a] -> Cycle a
fromList = foldr addMove mempty

-- Always check that the names are unique, because they are used to identify the
-- moves.
instance Semigroup (Cycle a) where
  l <> (Cycle xs) = foldr addMove l xs

instance Monoid (Cycle a) where
  mempty = Cycle []

-- | Change 'Move's in 'Cycle'.
mapCycle :: (Move a -> Move a) -> Cycle a -> Cycle a
mapCycle f = Cycle . map f . fromCycle

left :: Int -> String -> String
left n s = take n s ++ replicate (n - l) ' '
  where l = length s

right :: Int -> String -> String
right n s = replicate (n - l) ' ' ++ take n s
  where l = length s

renderRow :: String -> String -> String -> String -> String
renderRow name weight acceptRatio tuneParam = nB ++ wB ++ rB ++ tB
  where nB = left  25 name
        wB = right 8 weight
        rB = right 20 acceptRatio
        tB = right 20 tuneParam

moveHeader :: String
moveHeader = renderRow "Move name" "Weight" "Acceptance ratio" "Tuning parameter"

summarizeMove :: Move a -> Maybe Double -> String
summarizeMove m r = renderRow name weight acceptRatio tuneParamStr
  where name         = mvName m
        weight       = show $ mvWeight m
        acceptRatio  = maybe "" (printf "%.3f") r
        tuneParamStr = maybe "" (printf "%.3f") (tParam <$> mvTune m)

-- | Summarize the 'Move's in the 'Cycle'. Also report acceptance ratios for the
-- given number of last cycles.
summarizeCycle :: Cycle a -> String
summarizeCycle c = unlines $
  [
    replicate (length moveHeader) '─'
  , moveHeader
  , replicate (length moveHeader) '─'
  ] ++
  [ summarizeMove m Nothing | m <- mvs ] ++
  [ replicate (length moveHeader) '─'
  , show mpc ++ " move(s) per cycle." ]
  where mvs   = fromCycle c
        mpc   = sum $ map mvWeight mvs

-- | See 'summarizeCycle'. Also report acceptance ratios for the given number of
-- last cycles.
summarizeCycleA :: Int -> Acceptance (Move a) -> Cycle a -> String
summarizeCycleA n a c = unlines $
  [
    replicate (length moveHeader) '─'
  , moveHeader
  , replicate (length moveHeader) '─'
  ] ++
  [ summarizeMove m (ars M.!? m) | m <- mvs ] ++
  [ replicate (length moveHeader) '─'
  , show mpc ++ " move(s) per cycle." ++ arStr ]
  where mvs   = fromCycle c
        mpc   = sum $ map mvWeight mvs
        arStr = " Acceptance ratio(s) calculated over " ++ show n ++ " iterations."
        ars   = acceptanceRatios n a

-- | For each key @k@, store the list of accepted (True) and rejected (False)
-- proposals. For reasons of efficiency, the lists are stored in reverse order;
-- latest first.
newtype Acceptance k = Acceptance { fromAcceptance :: Map k [Bool] }

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
empty :: Ord k => [k] -> Acceptance k
empty ks = Acceptance $ M.fromList [ (k, []) | k <- ks ]

-- | For key @k@, prepend an accepted (True) or rejected (False) proposal.
prependA :: (Ord k, Show k) => k -> Bool -> Acceptance k -> Acceptance k
-- XXX: Unsafe; faster.
prependA k v (Acceptance m) = Acceptance $ M.adjust (v:) k m
-- -- XXX: Safe; slower.
-- prependA k v (Acceptance m) | k `M.member` m = Acceptance $ M.adjust (v:) k m
--                             | otherwise = error msg
--   where msg = "prependA: Can not add acceptance value for key: " <> show k <> "."
{-# INLINEABLE prependA #-}

-- | Reset acceptance storage.
resetA :: Acceptance k -> Acceptance k
resetA = Acceptance . M.map (const []) . fromAcceptance

ratio :: [Bool] -> Double
ratio xs = fromIntegral (length ts) / fromIntegral (length xs)
  where ts = filter (==True) xs

-- | Acceptance ratios averaged over the given number of last iterations. If
-- less than @n@ iterations are available, only those are used.
acceptanceRatios :: Int -> Acceptance k -> Map k Double
acceptanceRatios n (Acceptance m)= M.map (ratio . take n) m
