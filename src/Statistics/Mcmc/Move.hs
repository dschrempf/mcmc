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

-- TODO: Report moves; report tuning parameter and acceptance ratio.

module Statistics.Mcmc.Move
  (
    -- * Move
    Move (..)
  , MoveSimple (..)
  , Tuner
  , tuner
  , tune
  , autotune
    -- * Cycle
  , Cycle (fromCycle)
  , addMove
  , fromList
  , moves
  , mapCycle
  ) where

import Data.Bifunctor
import Data.Function
import Numeric.Log
import System.Random.MWC

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
data Tuner a = Tuner Double (Double -> MoveSimple a)

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
-- In detail, a 'Cycle' is a list of tuples @(move, weight)@. The weight
-- determines how often a 'Move' is executed per 'Cycle'.
--
-- 'Move's are replicated according to their weights and executed in random
-- order. That is, they are not executed in the order they appear in the
-- 'Cycle'.
newtype Cycle a = Cycle { fromCycle :: [(Move a, Int)] }
-- XXX: A classical list is used, and not 'Data.List.NonEmpty' (keep things
-- simple).

-- | Add a 'Move' with given weight to the 'Cycle'. The name of the added 'Move'
-- must be unique.
addMove :: Move a -> Int -> Cycle a -> Cycle a
addMove m w c | m `notElem` moves c = Cycle $ (m, w) : fromCycle c
              | otherwise = error msg
  where msg = "addMove: Move " <> mvName m <> " already exists in cycle."

-- | Create a 'Cycle' from a list of 'Move's with associated weights.
fromList :: [(Move a, Int)] -> Cycle a
fromList = foldr (uncurry addMove) mempty

-- Always check that the names are unique, because they are used to identify the
-- moves.
instance Semigroup (Cycle a) where
  l <> (Cycle xs) = foldr (uncurry addMove) l xs

instance Monoid (Cycle a) where
  mempty = Cycle []

-- | Extract 'Move's from 'Cycle'.
moves :: Cycle a -> [Move a]
moves = map fst . fromCycle

-- | Change 'Move's in 'Cycle'.
mapCycle :: (Move a -> Move a) -> Cycle a -> Cycle a
mapCycle f = Cycle . map (first f) . fromCycle
