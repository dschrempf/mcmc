{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc
-- Description :  Markov chain Monte Carlo samplers, batteries included
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 18:01:15 2020.
--
-- For an introduction to Markov chain Monte Carlo (MCMC) samplers and update
-- mechanisms using the Metropolis-Hastings-Green algorithm, please see Geyer,
-- C. J., 2011; Introduction to Markov Chain Monte Carlo. In Handbook of Markov
-- Chain Monte Carlo (pp. 45), Chapman \& Hall/CRC.
--
-- __The import of this module alone should cover most use cases.__
--
-- An MCMC sampler can be run with 'mcmc', for example using the
-- Metropolis-Hastings-Green algorithm 'mhg'.
--
-- The following example infers the mean deviation of a normally distributed
-- variable. For more involved inferences, please see
-- [mcmc-examples](https://github.com/dschrempf/mcmc/tree/master/mcmc-examples).
--
--
-- @
-- import Control.Monad
-- import Mcmc
-- import System.Random.MWC
--
-- trueMean, trueStdDev :: Double
-- trueMean = 5
-- trueStdDev = 4
--
-- lh :: LikelihoodFunction Double
-- lh = normal trueMean trueStdDev
--
-- cc :: Cycle Double
-- cc = cycleFromList [slideSymmetric 1.0 (PName "Medium") (PWeight 1) Tune]
--
-- mons :: [MonitorParameter Double]
-- mons = [monitorDouble "mu"]
--
-- monStd :: MonitorStdOut Double
-- monStd = monitorStdOut mons 200
--
-- mon :: Monitor Double
-- mon = Monitor monStd [] []
--
-- runMcmc :: GenIO -> IO (MHG Double)
-- runMcmc g = do
--   let s = Settings
--             (AnalysisName \"Normal\")
--             (BurnInWithAutoTuning 2000 200)
--             20000 Overwrite Sequential
--             NoSave Quiet
--       a = mhg noPrior lh cc mon 0 g
--   mcmc s a
-- @
module Mcmc
  ( -- * Proposals

    -- | A 'Proposal' is an instruction about how to advance a given Markov
    -- chain so that it possibly reaches a new state. That is, 'Proposal's
    -- specify how the chain traverses the state space. As far as this MCMC
    -- library is concerned, 'Proposal's are considered to be /elementary
    -- updates/ in that they cannot be decomposed into smaller updates.
    --
    -- 'Proposal's can be combined to form composite updates, a technique often
    -- referred to as /composition/. On the other hand, /mixing/ (used in the
    -- sense of mixture models) is the random choice of a 'Proposal' (or a
    -- composition of 'Proposal's) from a given set.
    --
    -- The __composition__ and __mixture__ of 'Proposal's allows specification
    -- of nearly all MCMC algorithms involving a single chain (i.e., population
    -- methods such as particle filters are excluded). In particular, Gibbs
    -- samplers of all sorts can be specified using this procedure. For
    -- reference, please see the short [encyclopedia of MCMC
    -- methods](https://dschrempf.github.io/coding/2020-11-12-encyclopedia-of-markov-chain-monte-carlo-methods/).
    --
    -- This library enables composition and mixture of 'Proposal's via the 'Cycle'
    -- data type. Essentially, a 'Cycle' is a set of 'Proposal's. The chain advances
    -- after the completion of each 'Cycle', which is called an __iteration__,
    -- and the iteration counter is increased by one.
    --
    -- The 'Proposal's in a 'Cycle' can be executed in the given order or in a
    -- random sequence which allows, for example, specification of a fixed scan
    -- Gibbs sampler, or a random sequence scan Gibbs sampler, respectively. See
    -- 'Order'.
    --
    -- Note that it is of utter importance that the given 'Cycle' enables
    -- traversal of the complete state space. Otherwise, the Markov chain will
    -- not converge to the correct stationary posterior distribution.
    --
    -- Proposals are named according to what they do, i.e., how they change the
    -- state of a Markov chain, and not according to the intrinsically used
    -- probability distributions. For example, 'slideSymmetric' is a sliding
    -- proposal. Under the hood, it uses the normal distribution with mean zero and
    -- given variance. The sampled variate is added to the current value of the
    -- variable (hence, the name slide). The same nomenclature is used by
    -- RevBayes [1]. The probability distributions and intrinsic properties of a
    -- specific proposal are specified in detail in the documentation.
    --
    -- The other method, which is used intrinsically, is more systematic, but
    -- also a little bit more complicated: we separate between the proposal
    -- distribution and how the state is affected. And here, I am referring to
    -- the operator (addition, multiplication, any other binary operator). For
    -- example, the sliding proposal with mean @m@, standard deviation @s@, and
    -- tuning parameter @t@ is implemented as
    --
    -- @
    -- slideSimple :: Double -> Double -> Double -> ProposalSimple Double
    -- slideSimple m s t =
    --   genericContinuous (normalDistr m (s * t)) (+) (Just negate) Nothing
    -- @
    --
    -- This specification is more involved. Especially since we need to know the
    -- probability of jumping back, and so we need to know the inverse operator
    -- 'negate'. However, it also allows specification of new proposals with
    -- great ease.
    --
    -- [1] Höhna, S., Landis, M. J., Heath, T. A., Boussau, B., Lartillot, N., Moore,
    -- B. R., Huelsenbeck, J. P., …, Revbayes: bayesian phylogenetic inference using
    -- graphical models and an interactive model-specification language, Systematic
    -- Biology, 65(4), 726–736 (2016). http://dx.doi.org/10.1093/sysbio/syw021
    PName (..),
    PWeight (..),
    Proposal,
    (@~),
    Tune (..),
    scale,
    scaleUnbiased,
    scaleContrarily,
    scaleBactrian,
    slide,
    slideSymmetric,
    slideUniformSymmetric,
    slideContrarily,
    slideBactrian,
    module Mcmc.Proposal.Simplex,
    Cycle,
    cycleFromList,
    Order (..),
    setOrder,

    -- * Settings
    Settings (..),
    -- ** Data types
    AnalysisName (..),
    BurnInSpecification (..),
    NIterations (..),
    ExecutionMode (..),
    ParallelizationMode (..),
    SaveMode (..),
    Verbosity (..),

    -- * Monitor

    -- | A 'Monitor' describes which part of the Markov chain should be logged
    -- and where. There are three different types:
    --
    -- ['MonitorStdOut'] Log to standard output.
    --
    -- ['MonitorFile'] Log to a file.
    --
    -- ['MonitorBatch'] Log summary statistics such as the mean of the last
    -- states to a file.
    Monitor (Monitor),
    MonitorStdOut,
    monitorStdOut,
    MonitorFile,
    monitorFile,
    MonitorBatch,
    monitorBatch,
    module Mcmc.Monitor.Parameter,
    module Mcmc.Monitor.ParameterBatch,

    -- * Prior distributions

    -- | Convenience functions for computing priors.
    module Mcmc.Prior,

    -- * Run and continue MCMC samplers
    mcmc,
    mcmcContinue,

    -- * Algorithms
    -- ** Metropolis-Hastings-Green algorithm
    MHG,
    mhg,
    -- ** Metropolis-coupled Markov chain Monte Carlo algorithm
    NChains (..),
    SwapPeriod (..),
    NSwaps (..),
    MC3Settings (..),
    MC3,
    mc3,

    -- * Save and load
    settingsLoad,
    mhgSave,
    mhgLoad,
    mc3Save,
    mc3Load,

    -- * Useful type synonyms
    PriorFunction,
    noPrior,
    LikelihoodFunction,
    noLikelihood,
  )
where

import Mcmc.Algorithm.MC3
import Mcmc.Algorithm.Metropolis
import Mcmc.Chain.Chain
import Mcmc.Mcmc
import Mcmc.Monitor
import Mcmc.Monitor.Parameter
import Mcmc.Monitor.ParameterBatch
import Mcmc.Prior
import Mcmc.Proposal
import Mcmc.Proposal.Bactrian
import Mcmc.Proposal.Scale
import Mcmc.Proposal.Simplex
import Mcmc.Proposal.Slide
import Mcmc.Settings
