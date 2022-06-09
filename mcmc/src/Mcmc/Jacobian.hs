-- |
-- Module      :  Mcmc.Jacobian
-- Description :  Types and convenience functions for Jacobians
-- Copyright   :  2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Tue May 31 09:37:54 2022.
module Mcmc.Jacobian
  ( Jacobian,
    JacobianG,
    JacobianFunction,
    JacobianFunctionG,
  )
where

import Numeric.Log

-- | Absolute value of the determinant of the Jacobian matrix.
type Jacobian = Log Double

-- | Generalized Jacobian.
type JacobianG a = Log a

-- | Function calculating the 'Jacobian'.
type JacobianFunction a = JacobianFunctionG a Double

-- | Function calculating the 'Jacobian'.
type JacobianFunctionG a b = a -> JacobianG b
