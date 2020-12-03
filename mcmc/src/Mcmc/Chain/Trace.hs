-- |
-- Module      :  Mcmc.Chain.Trace
-- Description :  History of a Markov chain
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 09:11:25 2020.
module Mcmc.Chain.Trace
  ( Trace,
    replicateT,
    pushT,
    headT,
    takeT,
    freezeT,
    thawT,
  )
where

import Control.Monad.Primitive
import qualified Data.Stack.Circular as C
import qualified Data.Vector as VB
import Mcmc.Chain.Link

-- | A 'Trace' is a mutable circular stack that passes through a list of states
-- with associated likelihoods which are called 'Link's.
newtype Trace a = Trace {fromTrace :: C.MStack VB.Vector RealWorld (Link a)}

-- | The empty trace.
replicateT :: Int -> Link a -> IO (Trace a)
replicateT n l = Trace <$> C.replicate n l

-- | Prepend an 'Link' to a 'Trace'.
pushT :: Link a -> Trace a -> IO (Trace a)
pushT x t = do
  s' <- C.push x (fromTrace t)
  return $ Trace s'
{-# INLINEABLE pushT #-}

-- | Get the most recent links of the trace.
--
-- O(1).
headT :: Trace a -> IO (Link a)
headT = C.get . fromTrace
{-# INLINEABLE headT #-}

-- | Get the k most recent links of the trace.
--
-- O(k).
takeT :: Int -> Trace a -> IO (VB.Vector (Link a))
takeT k = C.take k . fromTrace

-- | Freeze the mutable trace for storage.
freezeT :: Trace a -> IO (C.Stack VB.Vector (Link a))
freezeT = C.freeze . fromTrace

-- | Thaw a circular stack.
thawT :: C.Stack VB.Vector (Link a) -> IO (Trace a)
thawT t = Trace <$> C.thaw t
