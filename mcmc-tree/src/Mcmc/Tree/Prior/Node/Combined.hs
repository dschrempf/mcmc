-- |
-- Module      :  Mcmc.Tree.Prior.Node.Combined
-- Description :  Combined calibrations and constraints
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Fri Jun 25 11:23:24 2021.
module Mcmc.Tree.Prior.Node.Combined
  ( calibrateAndConstrain,
  )
where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import ELynx.Tree
import Mcmc.Chain.Chain
import Mcmc.Statistics.Types
import Mcmc.Tree.Prior.Node.Calibration
import Mcmc.Tree.Prior.Node.Constraint
import Mcmc.Tree.Types

-- Get the heights of all nodes and store them in a vector.
getAllHeights :: HasHeight a => Tree e a -> VU.Vector Height
getAllHeights = VU.fromList . map getHeight . labels

calibrateV :: StandardDeviation -> Calibration -> PriorFunction (VU.Vector Height)
calibrateV s c hs = calibrateSoftF s l h
  where
    l = calibrationInterval c
    i = calibrationNodeIndex c
    h = hs VU.! i

constrainV :: StandardDeviation -> Constraint -> PriorFunction (VU.Vector Height)
constrainV s k hs = constrainSoftF s (hY, hO)
  where
    iY = constraintYoungNodeIndex k
    hY = hs VU.! iY
    iO = constraintOldNodeIndex k
    hO = hs VU.! iO

-- | Calibrate and constrain nodes.
--
-- See 'calibrate', and 'constrain'.
--
-- This first extracts all node heights from the trees and then checks the
-- calibrations and constraints.
--
-- Use if you have many calibrations or constraints.
--
-- Do not use, if only a few calibrations and constraints have to be checked.
calibrateAndConstrain ::
  HasHeight a =>
  -- | Standard deviation of calibrations.
  StandardDeviation ->
  V.Vector Calibration ->
  -- | Height multiplier of tree for calibrations.
  Double ->
  -- | Standard deviation of constraints.
  StandardDeviation ->
  V.Vector Constraint ->
  PriorFunction (Tree e a)
calibrateAndConstrain sdC cs h sdK ks t
  | h <= 0 = error "calibrate: Height multiplier is zero or negative."
  | otherwise = V.product csPr * V.product ksPr
  where
    hs = getAllHeights t
    transform (Calibration n x i l) =
      let l' = if h == 1.0 then l else transformInterval (recip h) l
       in Calibration n x i l'
    csPr = V.map ((\c -> calibrateV sdC c hs) . transform) cs
    ksPr = V.map (\k -> constrainV sdK k hs) ks
