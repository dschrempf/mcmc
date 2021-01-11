-- |
-- Module      :  Mcmc.MarginalLikelihood
-- Description :  Calculate the marginal likelihood
-- Copyright   :  (c) Dominik Schrempf 2021
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Jan 11 16:34:18 2021.
module Mcmc.MarginalLikelihood
  ( NPoints (..),
    marginalLikelihood,
  )
where

-- TODO: Stepping stone sampling.
--
-- See Xie2010 and Fan2010.

import qualified Control.Monad.Parallel as P
import Data.Aeson
import qualified Data.Vector as VB
import Mcmc.Algorithm.MHG
import Mcmc.Chain.Chain
import Mcmc.Chain.Link
import Mcmc.Chain.Trace
import Mcmc.Internal.Random
import Mcmc.Mcmc
import Mcmc.Monitor
import Mcmc.Proposal
import Mcmc.Settings
import Numeric.Log
import System.Random.MWC
import Prelude hiding (cycle)

-- | The number of points used to approximate the path integral.
newtype NPoints = NPoints {fromNPoints :: Int}
  deriving (Eq, Read, Show)

-- See Figure 1 in Höhna, S., Landis, M. J., & Huelsenbeck, J. P., Parallel
-- power posterior analyses for fast computation of marginal likelihoods in
-- phylogenetics (2017). http://dx.doi.org/10.1101/104422.
getBetas :: NPoints -> [Log Double]
getBetas x = [f i ** (1.0 / 0.3) | i <- [0 .. k1]]
  where
    k = fromNPoints x
    k1 = pred k
    f j = Exp $ log $ fromIntegral j / fromIntegral k1

-- TODO: Check acceptance ratio and warn if low or high.
goToBeta ::
  ToJSON a =>
  Log Double ->
  BurnInSpecification ->
  Iterations ->
  LikelihoodFunction a ->
  MHG a ->
  IO (MHG a)
goToBeta b bi is lhf a = do
  mcmc ss a'
  where
    -- MCMC Settings.
    nm = AnalysisName "marginal-likelihood"
    ss = Settings nm bi is Fail Sequential NoSave Quiet
    -- Amend the likelihood function.
    lhf' = (** b) . lhf
    -- Amend the MHG algorithm.
    ch = fromMHG a
    l = link ch
    ch' =
      ch
        { -- Important: Update the likelihood using the new likelihood function.
          link = l {likelihood = lhf' $ state l},
          iteration = 0,
          start = 0,
          likelihoodFunction = lhf'
        }
    a' = MHG ch'

traverseBetas ::
  ToJSON a =>
  [Log Double] ->
  BurnInSpecification ->
  Iterations ->
  LikelihoodFunction a ->
  MHG a ->
  -- Posterior probabilities.
  IO [Log Double]
traverseBetas [] _ _ _ _ = return []
traverseBetas (b : bs) bi is lhf a = do
  -- Go to the next beta.
  a' <- goToBeta b bi is lhf a
  -- Extract the links.
  ls <- takeT n $ trace $ fromMHG a'
  -- Calculate the mean posterior probability.
  let mp = VB.sum (VB.map getPosterior ls) / fromIntegral (VB.length ls)
  -- Get the mean posterior probabilities of the other betas.
  mps <- traverseBetas bs bi is lhf a'
  return $ mp : mps
  where
    n = fromIterations is
    getPosterior l = prior l * likelihood l

-- TODO: Proper return value; marginal likelihood and confidence interval.

-- TODO: Proper output.

-- | Calculate the marginal likelihood using a path integral.
--
-- Also known as thermodynamic integration. In particular, /Annealing-Melting
-- Integration/ is used.
--
-- See Lartillot, N., & Philippe, H., Computing Bayes Factors Using
-- Thermodynamic Integration, Systematic Biology, 55(2), 195–207 (2006).
-- http://dx.doi.org/10.1080/10635150500433722
marginalLikelihood ::
  ToJSON a =>
  NPoints ->
  -- | Initial burn in at the starting point of the path.
  BurnInSpecification ->
  -- | Repetitive burn in at each point on the path.
  BurnInSpecification ->
  -- | The number of iterations performed at each point.
  Iterations ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  -- | Initial state.
  a ->
  GenIO ->
  -- TODO: Document return value.
  IO (Log Double, Log Double)
marginalLikelihood ps biI biR is prf lhf cc i0 g = do
  [g0, g1] <- splitGen 2 g
  [a0, a1] <-
    P.sequence
      [ mhg prf lhf cc mn trLen i0 g0,
        mhg prf lhf cc mn trLen i0 g1
      ]
  [mhg0, mhg1] <-
    P.sequence
      [ goToBeta 0.0 biI is lhf a0,
        goToBeta 1.0 biI is lhf a1
      ]
  [mps0, mps1] <- P.sequence
    [ traverseBetas bs0 biR is lhf mhg0,
      traverseBetas bs1 biR is lhf mhg1
    ]
  return (triangle mps0 bs0, triangle (reverse mps1) bs0)
  where
    trLen = TraceMinimum $ fromIterations is
    mn = noMonitor 1
    bs0 = getBetas ps
    bs1 = reverse bs0

triangle ::
  -- Y values.
  [Log Double] ->
  -- X values.
  [Log Double] ->
  -- Integral.
  Log Double
triangle (x0 : x1 : xs) (b0 : b1 : bs) = (x0 + x1) / (b1 - b0) + triangle (x1 : xs) (b1 : bs)
triangle _ _ = 0
