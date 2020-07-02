{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Mcmc
-- Description :  Markov chain Monte Carlo algorithms, batteries included
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue May  5 18:01:15 2020.
--
-- A short introduction to update mechanisms using the Metropolis-Hastings
-- algorithm (see Geyer, C. J., 2011; Introduction to Markov Chain Monte Carlo. In
-- Handbook of Markov Chain Monte Carlo (pp. 45), Chapman \& Hall/CRC).
--
-- For examples, please see
-- [mcmc-examples](https://github.com/dschrempf/mcmc/tree/master/mcmc-examples).
--
-- __The import of this module alone should cover most use cases.__
module Mcmc
  ( -- * Moves

    -- | A 'Move' is an instruction about how to advance a given Markov chain so
    -- that it possibly reaches a new state. That is, 'Move's specify how the
    -- chain traverses the state space. As far as this MCMC library is
    -- concerned, 'Move's are /elementary updates/ in that they cannot be
    -- decomposed into smaller updates.
    --
    -- 'Move's can be combined to form composite updates, a technique often
    -- referred to as /composition/. On the other hand, /mixing/ (used in the
    -- sense of mixture models) is the random choice of a 'Move' (or a
    -- composition of 'Move's) from a given set.
    --
    -- The __composition__ and __mixture__ of 'Move's allows specification of
    -- nearly all MCMC algorithms involving a single chain (i.e., population
    -- methods such as particle filters are excluded). In particular, Gibbs
    -- samplers of all sorts can be specified using this procedure.
    --
    -- This library enables composition and mixture of 'Move's via the 'Cycle'
    -- data type. Essentially, a 'Cycle' is a set of 'Move's. The chain advances
    -- after the completion of each 'Cycle', which is called an __iteration__,
    -- and the iteration counter is increased by one.
    --
    -- The 'Move's in a 'Cycle' can be executed in the given order or in a
    -- random sequence which allows, for example, specification of a fixed scan
    -- Gibbs sampler, or a random sequence scan Gibbs sampler, respectively. See
    -- 'Order'.
    --
    -- Note that it is of utter importance that the given 'Cycle' enables
    -- traversal of the complete state space. Otherwise, the Markov chain will
    -- not converge to the correct stationary posterior distribution.
    --
    -- Moves are named according to what they do, i.e., how they change the
    -- state of a Markov chain, and not according to the intrinsically used
    -- probability distributions. For example, 'slideSymmetric' is a sliding
    -- move. Under the hood, it uses the normal distribution with mean zero and
    -- given variance. The sampled variate is added to the current value of the
    -- variable (hence, the name slide). The same nomenclature is used by
    -- RevBayes [1]. The probability distributions and intrinsic properties of a
    -- specific move are specified in detail in the documentation.
    --
    -- The other method, which is used intrinsically, is more systematic, but
    -- also a little bit more complicated: we separate between the proposal
    -- distribution and how the state is affected. And here, I am not only
    -- referring to the accessor (i.e., the lens), but also to the operator
    -- (addition, multiplication, any other binary operator). For example, the
    -- sliding move (without tuning information) is implemented as
    --
    -- @
    -- slideSimple :: Lens' a Double -> Double -> Double -> Double -> MoveSimple a
    -- slideSimple l m s t = moveGenericContinuous l (normalDistr m (s * t)) (+) (-)
    -- @
    --
    -- This specification is more involved. Especially since we need to know the
    -- probability of jumping back, and so we need to know the inverse operator.
    -- However, it also allows specification of new moves with great ease.
    --
    -- [1] Höhna, S., Landis, M. J., Heath, T. A., Boussau, B., Lartillot, N., Moore,
    -- B. R., Huelsenbeck, J. P., …, Revbayes: bayesian phylogenetic inference using
    -- graphical models and an interactive model-specification language, Systematic
    -- Biology, 65(4), 726–736 (2016). http://dx.doi.org/10.1093/sysbio/syw021
    Move,
    scale,
    scaleUnbiased,
    slide,
    slideBactrian,
    slideSymmetric,
    slideUniform,
    Cycle,
    fromList,
    Order (..),
    setOrder,

    -- * Initialization

    -- | The 'Status' contains all information to run an MCMC chain. It is
    -- constructed using the function 'status'.
    --
    -- The 'Status' of a Markov chain includes information about current state
    -- ('Mcmc.Item.Item') and iteration, the history of the chain
    -- ('Mcmc.Trace.Trace'), the 'Acceptance' ratios, and the random number
    -- generator.
    --
    -- Further, the 'Status' includes auxiliary variables and functions such as
    -- the prior and likelihood functions, instructions to move around the state
    -- space (see above) and to monitor the MCMC run, as well as some auxiliary
    -- information.
    status,
    noSave,
    quiet,
    debug,

    -- * Monitor

    -- | A 'Monitor' describes which part of the Markov chain should be logged
    -- and where. There are three different types:
    -- - 'MonitorStdOut': Log to standard output.
    -- - 'MonitorFile': Log to a file.
    -- - 'MonitorBatch': Log summary statistics such as the mean of the last
    -- - states to a file.
    Monitor (Monitor),
    MonitorStdOut,
    monitorStdOut,
    MonitorFile,
    monitorFile,
    MonitorBatch,
    monitorBatch,
    module Mcmc.Monitor.Parameter,
    module Mcmc.Monitor.ParameterBatch,

    -- * Algorithms

    -- | At the moment, the library is tailored to the Metropolis-Hastings
    -- algorithm ('mh') since it covers most use cases. However, implementation
    -- of more algorithms is planned in the future.
    mh,
    mhContinue,

    -- * Misc
    pzero,
    loadStatus,
  )
where

import Mcmc.Metropolis
import Mcmc.Monitor
import Mcmc.Monitor.Parameter
import Mcmc.Monitor.ParameterBatch
import Mcmc.Move
import Mcmc.Move.Bactrian
import Mcmc.Move.Scale
import Mcmc.Move.Slide
import Mcmc.Save
import Mcmc.Status
import Numeric.Log

-- | Because we need a probability of zero for likelihoods of really bad moves.
pzero :: Fractional a => Log a
pzero = Exp $ - (1 / 0)
