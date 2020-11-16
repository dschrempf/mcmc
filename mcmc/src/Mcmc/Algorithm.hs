-- |
-- Module      :  Mcmc.Algorithm
-- Description :  MCMC algorithms
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Nov 16 14:37:11 2020.
module Mcmc.Algorithm
  ( Algorithm (..),
  )
where

import Numeric.Log

-- | TODO: REFACTOR. Documentation.
class Algorithm a where
  -- | Get current iteration.
  getIteration :: a -> Int

  -- | Perform one iteration.
  jump :: a -> a

  -- | Auto tune the proposals.
  autotune :: a -> a

  -- | Clean the state.
  clean :: a -> a

  -- TODO: REFACTOR. 'Save' and 'loadWith'.

  -- | Report prior and likelihood; useful for debugging.
  report :: a -> (Log Double, Log Double)

-- -- | Auto tune the 'Proposal's in the 'Cycle' of the chain. Reset acceptance counts.
-- -- See 'autotuneCycle'.
-- mcmcAutotune :: Mcmc a ()
-- mcmcAutotune = do
--   mcmcDebugB "Auto tune."
--   s <- get
--   let a = acceptance s
--       c = cycle s
--       c' = autotuneCycle a c
--   put $ s {cycle = c'}

-- -- | Clean the state.
-- mcmcClean :: Mcmc a ()
-- mcmcClean = do
--   s <- get
--   let cl = cleaner s
--       i = iteration s
--   case cl of
--     Just (Cleaner n f) | i `mod` n == 0 -> do
--       mcmcDebugB "Clean state."
--       let (Item st pr lh) = item s
--       mcmcDebugS $
--         "Old log prior and log likelihood: " ++ show (ln pr) ++ ", " ++ show (ln lh) ++ "."
--       let prF = priorF s
--           lhF = likelihoodF s
--           st' = f st
--           pr' = prF st'
--           lh' = lhF st'
--       mcmcDebugS $
--         "New log prior and log likelihood: " ++ show (ln pr') ++ ", " ++ show (ln lh') ++ "."
--       let dLogPr = abs $ ln pr - ln pr'
--           dLogLh = abs $ ln lh - ln lh'
--       when
--         (dLogPr > 0.01)
--         (mcmcWarnS $ "Log of old and new prior differ by " ++ show dLogPr ++ ".")
--       when
--         (dLogPr > 0.01)
--         (mcmcWarnS $ "Log of old and new likelihood differ by " ++ show dLogLh ++ ".")
--       put $ s {item = Item st' pr' lh'}
--     _ -> return ()
