{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Mcmc.Algorithm.Metropolis
-- Description :  Metropolis-Hastings-Green algorithm
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 20:11:30 2020.
--
-- Metropolis-Hastings algorithm.
module Mcmc.Algorithm.Metropolis
  ( MHG (..),
  )
where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.RWS.CPS
import Data.Aeson
import Data.Maybe
import Mcmc.Algorithm
import Mcmc.Chain.Chain
import Mcmc.Environment
import Mcmc.Mcmc
import Mcmc.Proposal
import Numeric.Log
import System.Random.MWC
import Prelude hiding (cycle)

newtype MHG a = MHG {fromMHG :: Chain a}

instance Algorithm MHG

-- The Metropolis-Hastings ratio.
--
-- 'Infinity' if fX is zero. In this case, the proposal is always accepted.
--
-- 'NaN' if (fY or q) and fX are zero. In this case, the proposal is always
-- rejected.

-- There is a discrepancy between authors saying that one should (a) always
-- accept the new state when the current posterior is zero (Chapter 4 of the
-- Handbook of Markov Chain Monte Carlo), or (b) almost surely reject the
-- proposal when either fY or q are zero (Chapter 1). Since I trust the author
-- of Chapter 1 (Charles Geyer) I choose to follow option (b).
mhgRatio :: Log Double -> Log Double -> Log Double -> Log Double -> Log Double
-- q = qYX / qXY * jXY; see 'ProposalSimple'.
-- j = Jacobian.
mhgRatio fX fY q j = fY / fX * q * j
{-# INLINE mhgRatio #-}

mhgPropose :: Proposal a -> Mcmc a ()
mhgPropose m = do
  let p = pSimple m
  c <- get
  let (Item x pX lX) = item c
      pF = priorF c
      lF = likelihoodF c
      a = acceptance c
      g = generator c
  -- 1. Sample new state.
  (!y, !q, !j) <- liftIO $ p x g
  -- 2. Calculate Metropolis-Hastings ratio.
  let !pY = pF y
      !lY = lF y
      !r = mhRatio (pX * lX) (pY * lY) q j
  -- 3. Accept or reject.
  if ln r >= 0.0
    then do
      let !a' = pushA m True a
      put $ c {item = Item y pY lY, acceptance = a'}
    else do
      b <- uniform g
      if b < exp (ln r)
        then do
          let !a' = pushA m True a
          put $ c {item = Item y pY lY, acceptance = a'}
        else do
          let !a' = pushA m False a
          put $ c {acceptance = pushA m False a'}

-- TODO: Splitmix. Split the generator here. See SaveSpec -> mhContinue.

-- Run one iterations; perform all proposals in a Cycle.
mhgIter :: ToJSON a => [Proposal a] -> Mcmc a ()
mhgIter ps = do
  mapM_ mhPropose ps
  s <- get
  let i = item s
      t = trace s
      n = iteration s
  put $ s {trace = pushT i t, iteration = succ n}

-- Run N iterations.
mhgNIter :: ToJSON a => Int -> Mcmc a ()
mhgNIter n = do
  c <- gets cycle
  g <- gets generator
  cycles <- liftIO $ getNIterations c n g
  forM_ cycles mhIter
