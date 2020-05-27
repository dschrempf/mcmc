{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Statistics.Mcmc
Description :  Markov chain Monte Carlo algorithms, batteries included
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 18:01:15 2020.

-}

module Statistics.Mcmc
  (
    -- * Moves
    --
    -- Moves are named according to what they do, i.e., how they change the state of a
    -- Markov chain, and not according to the intrinsically used probability
    -- distributions. For example, 'slideDouble' is a sliding move changing a 'Double'.
    -- Under the hood, it uses the normal distribution with a given mean and variance.
    -- The sampled variate is added to the current value of the variable (hence, the
    -- name slide). The same nomenclature is used by RevBayes [1]. The probability
    -- distributions and intrinsic properties of a specific move are specified in
    -- detail in the documentation.
    --
    -- The other method, which is used intrinsically, is more systematic, but also a
    -- little bit more complicated: we separate between the proposal distribution and
    -- how the state is affected. And here, I am not only referring to the accessor
    -- (i.e., the lens), but also to the operator (addition, multiplication, any other
    -- binary operator). For example, the sliding move is implemented as
    --
    -- @
    -- slide l n m s = moveGenericContinuous l n (normalDistr m s) (+) (-)
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
    module Statistics.Mcmc.Move
  , module Statistics.Mcmc.Move.Slide
  , module Statistics.Mcmc.Move.Scale
    -- * Initialization
  , module Statistics.Mcmc.Status
  , module Statistics.Mcmc.Monitor
    -- * Algorithms
  , module Statistics.Mcmc.Metropolis
  ) where

import Statistics.Mcmc.Metropolis
import Statistics.Mcmc.Monitor
import Statistics.Mcmc.Move
import Statistics.Mcmc.Move.Scale
import Statistics.Mcmc.Move.Slide
import Statistics.Mcmc.Status

-- TODO: Moves on simplices: SimplexElementScale (?).

-- TODO: Moves on tree branch lengths.
-- - Slide a node on the tree.
-- - Scale a tree.

-- TODO: Moves on tree topologies.
-- - NNI
-- - Narrow (what is this, see RevBayes)
-- - FNPR (dito)

-- TODO: Bactrian moves; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3845170/.
--
-- slideBactrian
--
-- scaleBactrian
