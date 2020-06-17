{-# LANGUAGE RankNTypes #-}

{- |
Module      :  Mcmc
Description :  Markov chain Monte Carlo algorithms, batteries included
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue May  5 18:01:15 2020.

-}

module Mcmc
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
    module Mcmc.Move
  , module Mcmc.Move.Slide
  , module Mcmc.Move.Scale
    -- * Initialization
  , module Mcmc.Status
  , module Mcmc.Monitor
  , module Mcmc.Monitor.Parameter
  , module Mcmc.Monitor.ParameterBatch
    -- * Algorithms
  , module Mcmc.Metropolis
  )
where

import           Mcmc.Metropolis
import           Mcmc.Monitor
import           Mcmc.Monitor.Parameter
import           Mcmc.Monitor.ParameterBatch
import           Mcmc.Move
import           Mcmc.Move.Scale
import           Mcmc.Move.Slide
import           Mcmc.Status
