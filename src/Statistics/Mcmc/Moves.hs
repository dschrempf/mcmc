{- |
Module      :  Statistics.Mcmc.Moves
Description :  A collection of predefined moves
Copyright   :  (c) Dominik Schrempf 2020
License     :  GPL-3

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Thu May 14 13:51:51 2020.

I decided to provide moves named according to what they do, i.e., how they
change the state of a Markov chain, and not according to the intrinsically used
probability distributions. For example, 'moveSlide mean variance' acts on a
'Double', and uses the normal distribution with given @mean@ and @variance@. The
sampled variate is added to the current value of the variable (hence, the name
slide). The same nomenclature is used by RevBayes [1]. The probability
distributions and intrinsic properties of a specific move are specified in
detail in the respective help.

The other method is more systematic, but also a little bit more complicated: we
differentiate between the proposal distribution and how the state is affected.
Here, I am not only referring to the accessor (i.e., the lens), but also to the
operator (addition, multiplication, any other binary operator). In this scheme,
'mvSlide' would correspond to something like @move (normalDistribution mean
variance) (+) ...@. This specification is more involved. Especially, we need to
know the probability of jumping back, and so we need to know the reverse of the
operator. However, it also allows specification of new moves with great ease.

TODO: Implement a Common module providing @move distribution operator inverse@
and rewrite mvSlide to use this operator (check speed with benchmark).

TODO: Scaling move (Gamma distribution with mean 1.0).

TODO: Moves on simplices: SimplexElementScale (?).

TODO: Moves on tree branch lengths.
- Slide a node on the tree.
- Scale a tree.

TODO: Moves on tree topologies.
- NNI
- Narrow (what is this, see RevBayes)
- FNPR (dito)

[1] Höhna, S., Landis, M. J., Heath, T. A., Boussau, B., Lartillot, N., Moore,
B. R., Huelsenbeck, J. P., …, Revbayes: bayesian phylogenetic inference using
graphical models and an interactive model-specification language, Systematic
Biology, 65(4), 726–736 (2016). http://dx.doi.org/10.1093/sysbio/syw021
-}

module Statistics.Mcmc.Moves
  ( module Statistics.Mcmc.Moves.Slide
  ) where

import Statistics.Mcmc.Moves.Slide
