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

import qualified Data.Vector as VB
import ELynx.Tree
import Mcmc.Prior
import Mcmc.Statistics.Types
import Mcmc.Tree.Prior.Node.Calibration
import Mcmc.Tree.Prior.Node.Constraint
import Mcmc.Tree.Types

-- Get the heights of all nodes and store them in a vector.
getAllHeights :: HeightTree a -> VB.Vector a
getAllHeights = VB.fromList . branches . fromHeightTree

calibrateV ::
  (RealFloat a) =>
  StandardDeviation a ->
  Calibration a ->
  PriorFunctionG (VB.Vector a) a
calibrateV s c hs = calibrateSoftF s l h
  where
    l = calibrationInterval c
    i = calibrationNodeIndex c
    h = hs VB.! i

constrainV ::
  (RealFloat a) =>
  StandardDeviation a ->
  Constraint ->
  PriorFunctionG (VB.Vector a) a
constrainV s k hs = constrainSoftF s (hY, hO)
  where
    iY = constraintYoungNodeIndex k
    hY = hs VB.! iY
    iO = constraintOldNodeIndex k
    hO = hs VB.! iO

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
  (RealFloat a) =>
  -- | Standard deviation of calibrations.
  StandardDeviation a ->
  VB.Vector (Calibration a) ->
  -- | Height multiplier of tree for calibrations.
  a ->
  -- | Standard deviation of constraints.
  StandardDeviation a ->
  VB.Vector Constraint ->
  PriorFunctionG (HeightTree a) a
calibrateAndConstrain sdC cs h sdK ks t
  | h <= 0 = error "calibrate: Height multiplier is zero or negative."
  | otherwise = VB.product csPr * VB.product ksPr
  where
    hs = getAllHeights t
    transform (Calibration n x i l) =
      let l' = if h == 1.0 then l else transformInterval (recip h) l
       in Calibration n x i l'
    csPr = VB.map ((\c -> calibrateV sdC c hs) . transform) cs
    ksPr = VB.map (\k -> constrainV sdK k hs) ks
