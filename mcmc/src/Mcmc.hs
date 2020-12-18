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
-- C. J., (2011), Introduction to Markov Chain Monte Carlo, In Handbook of
-- Markov Chain Monte Carlo (pp. 45), CRC press.
--
-- This library focusses on classical Markov chain Monte Carlo algorithms such
-- as the Metropolis-Hastings-Green [1] algorithm, or population methods
-- involving parallel chains such as the Metropolic-coupled Markov chain Monte
-- Carlo [2] algorithm. In particular, sequential Monte Carlo [3] algorithms
-- following a moving posterior distribution are not provided.
--
-- An MCMC sampler can be run with 'mcmc', for example using the
-- Metropolis-Hastings-Green algorithm 'mhg'.
--
-- Usually, it is best to start with an example:
--
-- - Basic inference of the [accuracy of an
--   archer](https://github.com/dschrempf/mcmc/tree/master/mcmc-examples/Archery/Archery.hs)
--
-- - [More involved
--   examples](https://github.com/dschrempf/mcmc/tree/master/mcmc-examples/Archery/Archery.hs)
--
-- __The import of this module alone should cover most use cases.__
--
-- @[1]@ Geyer, C. J. (2011), Introduction to markov chain monte carlo, In
-- Handbook of Markov Chain Monte Carlo (pp. 45), CRC press.
--
-- @[2]@ Geyer, C. J. (1991), Markov chain monte carlo maximum likelihood,
-- Computing Science and Statistics, Proceedings of the 23rd Symposium on the
-- Interface.
--
-- @[3]@ Sequential monte carlo methods in practice (2001), Editors: Arnaud
-- Doucet, Nando de Freitas, and Neil Gordon, Springer New York.
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
    -- This library enables composition and mixture of 'Proposal's via the
    -- 'Cycle' data type. Essentially, a 'Cycle' is a set of 'Proposal's. The
    -- chain advances after the completion of each 'Cycle', which is called an
    -- __iteration__, and the iteration counter is increased by one.
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
    -- proposal. Under the hood, it uses the normal distribution with mean zero
    -- and given variance. The sampled variate is added to the current value of
    -- the variable (hence, the name slide). The same nomenclature is used by
    -- RevBayes [4]. The probability distributions and intrinsic properties of a
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
    -- @[4]@ Höhna, S., Landis, M. J., Heath, T. A., Boussau, B., Lartillot, N.,
    -- Moore, B. R., Huelsenbeck, J. P., …, Revbayes: bayesian phylogenetic
    -- inference using graphical models and an interactive model-specification
    -- language, Systematic Biology, 65(4), 726–736 (2016).
    -- http://dx.doi.org/10.1093/sysbio/syw021
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
    module Mcmc.Settings,

    -- * Monitors

    -- | A 'Monitor' describes which part of the Markov chain should be logged
    -- and where. Monitor files can directly be loaded into established MCMC
    -- analysis tools working with tab separated tables (for example,
    -- [Tracer](http://tree.bio.ed.ac.uk/software/tracer/)).
    --
    -- There are three different 'Monitor' types:
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

    -- * MCMC samplers
    mcmc,
    mcmcContinue,

    -- | See also 'settingsLoad', 'mhgLoad', and 'mc3Load'.

    -- * Algorithms
    module Mcmc.Algorithm.Metropolis,
    module Mcmc.Algorithm.MC3,


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
