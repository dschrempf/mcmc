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
    MLSettings (..),
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
import Text.Printf
import Prelude hiding (cycle)

-- Reciprocal temperature value traversed along the path integral.
type Point = Log Double

-- List of reciprocal temperature values traversed along the path integral.
type Points = [Point]

-- | The number of points used to approximate the path integral.
newtype NPoints = NPoints {fromNPoints :: Int}
  deriving (Eq, Read, Show)

data MLSettings = MLSettings
  { mlAnalysisName :: AnalysisName,
    mlNPoints :: NPoints,
    -- | Initial burn in at the starting point of the path.
    mlInitialBurnIn :: BurnInSpecification,
    -- | Repetitive burn in at each point on the path.
    mlPointBurnIn :: BurnInSpecification,
    -- | The number of iterations performed at each point.
    mlIterations :: Iterations
  }
  deriving (Eq, Read, Show)

-- Distribute the points according to a skewed beta distribution (more points at
-- low values). It is inconvenient that the reciprocal temperatures are denoted
-- as beta, and we also use the beta distribution :). Don't mix them up!
--
-- See discussion in Xie, W., Lewis, P. O., Fan, Y., Kuo, L., & Chen, M.,
-- Improving marginal likelihood estimation for bayesian phylogenetic model
-- selection, Systematic Biology, 60(2), 150–160 (2010).
-- http://dx.doi.org/10.1093/sysbio/syq085
--
-- Or Figure 1 in Höhna, S., Landis, M. J., & Huelsenbeck, J. P., Parallel power
-- posterior analyses for fast computation of marginal likelihoods in
-- phylogenetics (2017). http://dx.doi.org/10.1101/104422
getPoints :: NPoints -> Points
getPoints x = [f i ** (1.0 / 0.3) | i <- [0 .. k1]]
  where
    k = fromNPoints x
    k1 = pred k
    f j = Exp $ log $ fromIntegral j / fromIntegral k1

-- TODO: Check acceptance ratio and warn if low or high.
sampleAtPoint ::
  ToJSON a =>
  Point ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  IO (MHG a)
sampleAtPoint b ss lhf a = do
  mcmc ss' a'
  where
    -- For debugging set a proper analysis name.
    nm = sAnalysisName ss
    getName :: Point -> AnalysisName
    getName x = nm <> AnalysisName (printf "%.5f" $ exp $ ln x)
    ss' = ss { sAnalysisName = getName b}
    -- Amend the likelihood function. Don't calculate the likelihood when beta
    -- is 0.0.
    lhf' = if b == 0.0 then const 1.0 else (** b) . lhf
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

traversePoints ::
  ToJSON a =>
  Points ->
  Settings ->
  LikelihoodFunction a ->
  MHG a ->
  -- Posterior probabilities.
  IO [Log Double]
traversePoints [] _ _ _ = return []
traversePoints (b : bs) ss lhf a = do
  -- Go to the next beta.
  a' <- sampleAtPoint b ss lhf a
  -- Extract the links.
  ls <- takeT n $ trace $ fromMHG a'
  -- Calculate the mean posterior probability.
  let mp = VB.sum (VB.map getPosterior ls) / fromIntegral (VB.length ls)
  -- Get the mean posterior probabilities of the other betas.
  mps <- traversePoints bs ss lhf a'
  return $ mp : mps
  where
    n = fromIterations $ sIterations ss
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
  MLSettings ->
  PriorFunction a ->
  LikelihoodFunction a ->
  Cycle a ->
  -- | Initial state.
  a ->
  GenIO ->
  -- TODO: Document return value.
  IO (Log Double, Log Double)
marginalLikelihood env prf lhf cc i0 g = do
  [g0, g1] <- splitGen 2 g
  [a0, a1] <-
    P.sequence
      [ mhg prf lhf cc mn trLen i0 g0,
        mhg prf lhf cc mn trLen i0 g1
      ]
  [mhg0, mhg1] <-
    P.sequence
      [ sampleAtPoint 0.0 ssI lhf a0,
        sampleAtPoint 1.0 ssI lhf a1
      ]
  [mps0, mps1] <-
    P.sequence
      [ traversePoints bs0 ssP lhf mhg0,
        traversePoints bs1 ssP lhf mhg1
      ]
  return (triangle mps0 bs0, triangle (reverse mps1) bs0)
  where
    nm = mlAnalysisName env
    is = mlIterations env
    biI = mlInitialBurnIn env
    ssI = Settings nm biI is Fail Sequential NoSave Quiet
    biP = mlPointBurnIn env
    ssP = Settings nm biP is Fail Sequential NoSave Quiet
    trLen = TraceMinimum $ fromIterations is
    mn = noMonitor 1
    bs0 = getPoints $ mlNPoints env
    bs1 = reverse bs0

triangle ::
  -- Y values.
  [Log Double] ->
  -- X values.
  Points ->
  -- Integral.
  Log Double
triangle (x0 : x1 : xs) (b0 : b1 : bs) = (x0 + x1) / (b1 - b0) + triangle (x1 : xs) (b1 : bs)
triangle _ _ = 0
