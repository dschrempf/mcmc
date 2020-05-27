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

-}

module Statistics.Mcmc.Move
  (
    -- * Move
    MoveSimple (..)
  , Tuner (tParam, tFunct)
  , tuner
  , Move (..)
  , tune
  , autotune
    -- * Cycle
  , Cycle (fromCycle)
  , addMove
  , fromList
  , moves
  , mapC
  ) where

import Data.Bifunctor
import Data.Function
import Numeric.Log
import System.Random.MWC

-- | Simple move without tuning information.
--
-- We need to know the probability density of jumping forth, but also the
-- probability density of jumping back. They are needed to calculate the
-- Metropolis-Hastings ratio.
--
-- One could also use a different type for 'mvSample', so that 'mvLogDensity'
-- can be avoided. In detail,
--
-- @
--   mvSample :: a -> GenIO -> IO (a, Log Double, Log, Double)
-- @
--
-- where the log densities describe the probability of going there and back.
-- However, we may need more information about the move for other MCMC samplers
-- different from Metropolis-Hastings.
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
    tParam :: Double                 -- ^ Tuning parameter.
  , tFunct :: Double -> MoveSimple a -- ^ Tuning function; get a simple move for
                                     -- a specific tuning parameter.
  }

-- | Create a 'Tuner'.
tuner :: (Double -> MoveSimple a) -> Tuner a
tuner = Tuner 1.0

-- | A 'Move' is an instruction about how the Markov chain will traverse the
-- state space @a@. Essentially, it is a probability density conditioned on the
-- current state.
--
-- A 'Move' contains information about whether it is tuneable or not.
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
-- In detail, a 'Cycle' is list of tuples @(move, weight)@. The weight
-- determines how often a 'Move' is executed per 'Cycle'.
--
-- 'Move's are replicated according to their weights and executed in random
-- order. That is, they are not executed in the order they appear in the
-- 'Cycle'.
--
-- A classical list is used, and not 'Data.List.NonEmpty' (keep things simple).
newtype Cycle a = Cycle { fromCycle :: [(Move a, Int)] }

-- | Add a 'Move' with weight @w@ to the 'Cycle'. The name of the added 'Move'
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
mapC :: (Move a -> Move a) -> Cycle a -> Cycle a
mapC f = Cycle . map (first f) . fromCycle
