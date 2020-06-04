{- |
Module      :  Mcmc.Statistics.Statistics
Description :  Compute statistics about MCMC runs
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Wed Jun  3 18:38:18 2020.

-}

module Mcmc.Statistics.Statistics
  (
  ) where

-- TODO: This model is supposed to provide function computing summary statistics
-- of an MCMC run.
--
-- For example, I can think of the ESS, lag-k autocovariance and autocorrelation
-- (p 9), sample mean and sample variance, asymptotic variance (p 8).
--
-- Page number from Geyer, C. J., Introduction to markov chain monte carlo, In
-- Handbook of Markov Chain Monte Carlo (pp. 45) (2011). Chapman & Hall/CRC.
