{- |
Module      :  Mcmc.Save
Description :  Save the state of a Markov chain
Copyright   :  (c) Dominik Schrempf, 2020
License     :  GPL-3.0-or-later

Maintainer  :  dominik.schrempf@gmail.com
Stability   :  unstable
Portability :  portable

Creation date: Tue Jun 16 10:18:54 2020.

-}

-- TODO: Continue an MCMC run. It is easy to save and restore the current state
-- and likelihood (or the trace), but it is not feasible to store all the moves
-- and so on, so they have to be provided again when continuing a run.

-- TODO: upon continuation: recompute and check the posterior for the last state
-- because the posterior function may have changed. Of course, we cannot test
-- for the same function, but having the same posterior at the last state is
-- already a good indicator.


module Mcmc.Save
  (
  ) where

-- import           Data.Aeson

-- data Save = Save (Item a) Int (Trace a) ACCEPTANCE GEN

-- TODO THIS DOES NOT WORK! AN INTERMEDIATE DATA STRUCTURE IS NECESSARY. THIS
-- AGAIN, RAISES THE QUESTION; SHOULD I COMBINE THE TWO AGAIN??
--
--
--
-- instance (ToJSON a) => ToJSON (Status a) where
--   toJSON (Status s i t a g) = object
--     [
--       "s" .= s
--     , "i" .= i
--     , "t" .= t
--     , "a" .= a
--     , "g" .= g -- TODO
--     ]
--   toEncoding = undefined

-- instance (FromJSON a) => FromJSON (Status a) where
--   parseJSON = undefined
