-- |
-- Module      :  Mcmc.Chain.Trace
-- Description :  History of a Markov chain
-- Copyright   :  2021 Dominik Schrempf
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
    fromVectorT,
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

-- NOTE: We directly refer to the 'PrimSate' 'RealWorld' because, otherwise, we
-- have two type variables in 'Chain m a'.

-- | A 'Trace' is a mutable circular stack that passes through a list of states
-- with associated priors and likelihoods called 'Link's.
newtype Trace a = Trace {fromTrace :: C.MStack VB.Vector RealWorld (Link a)}

-- | Initialize a trace of given length by replicating the same value.
--
-- Be careful not to compute summary statistics before pushing enough values.
--
-- Call 'error' if the maximum size is zero or negative.
replicateT :: Int -> Link a -> IO (Trace a)
replicateT n l = Trace <$> C.replicate n l

-- | Create a trace from a vector. The length is determined by the vector.
--
-- Call 'error' if the vector is empty.
fromVectorT :: VB.Vector (Link a) -> IO (Trace a)
fromVectorT xs = Trace <$> C.fromVector xs

-- | Push a 'Link' on the 'Trace'.
pushT :: Link a -> Trace a -> IO (Trace a)
pushT x t = do
  s' <- C.push x (fromTrace t)
  return $ Trace s'
{-# INLINEABLE pushT #-}

-- | Get the most recent link of the trace (see 'C.get').
headT :: Trace a -> IO (Link a)
headT = C.get . fromTrace
{-# INLINEABLE headT #-}

-- | Get the k most recent links of the trace (see 'C.take').
takeT :: Int -> Trace a -> IO (VB.Vector (Link a))
takeT k = C.take k . fromTrace

-- | Freeze the mutable trace for storage (see 'C.freeze').
freezeT :: Trace a -> IO (C.Stack VB.Vector (Link a))
freezeT = C.freeze . fromTrace

-- | Thaw a circular stack (see 'C.thaw').
thawT :: C.Stack VB.Vector (Link a) -> IO (Trace a)
thawT t = Trace <$> C.thaw t
