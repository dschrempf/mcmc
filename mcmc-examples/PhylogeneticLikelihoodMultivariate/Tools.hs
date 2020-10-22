{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Tools
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Thu Oct 22 16:53:59 2020.
module Tools
  ( getBranches,
    sumFirstTwo,
    prepare,
  )
where

import Control.Monad
import Data.Aeson
import Data.Bifunctor
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List
import Data.Maybe
import qualified Data.Set as S
import qualified Data.Vector.Storable as V
import qualified ELynx.Topology.Rooted as T
import ELynx.Tree
import Mcmc.Tree
import qualified Numeric.LinearAlgebra as L

-- | Get all branches of a rooted tree. Store the branches in a vector such that
-- the two branches leading to the root are the first two entries of the vector.
-- Ignore the root branch.
getBranches :: Tree Double a -> V.Vector Double
getBranches (Node _ _ [l, r]) = V.fromList $ head ls : head rs : tail ls ++ tail rs
  where
    ls = branches l
    rs = branches r
getBranches _ = error "getBranches: Root node is not bifurcating."

-- | Sum the first two elements of a vector. Needed to merge the two branches
-- leading to the root.
sumFirstTwo :: V.Vector Double -> V.Vector Double
sumFirstTwo v = (v V.! 0 + v V.! 1) `V.cons` V.drop 2 v

-- Get the posterior matrix of branch lengths of rooted trees. Merge the two
-- branch lengths leading to the root.
getPosteriorMatrixRooted :: [Tree Double a] -> L.Matrix Double
getPosteriorMatrixRooted = L.fromRows . map (sumFirstTwo . getBranches)

-- Get the posterior matrix of branch lengths of unrooted trees (trees with
-- multifurcating root nodes). Before midpoint rooting, the mean branch lengths
-- of the unrooted trees have to be determined.
getPosteriorMatrix :: [Tree Double a] -> L.Matrix Double
getPosteriorMatrix = L.fromRows . map (V.fromList . branches)

-- | Read trees and extract branch lengths.
--
-- @prepare intrees bnAnalysis@
prepare :: FilePath -> String -> IO ()
prepare fnInTrees bnAnalysis = do
  putStrLn "Read trees; skip a burn in of 1000 trees."
  trs <- drop 1000 <$> someTrees fnInTrees

  putStrLn "Check if trees have the same topology."
  let l = length $ nub $ map T.fromLabeledTree trs
  unless (l == 1) (error "Trees have different topologies.")

  putStrLn "Calculate mean branch lengths."
  let pm = getPosteriorMatrix trs
      (means, _) = L.meanCov pm
  putStrLn "The mean branch lengths are:"
  print means
  let meanTree =
        fromMaybe (error "Could not label tree with mean branch lengths") $
          setBranches (map Length $ V.toList means) (head trs)
      lvs = leaves meanTree
      trOutgroup = either error id $ outgroup (S.singleton $ head lvs) "root" meanTree
      tr = either error id $ midpoint trOutgroup
  putStrLn "The tree with mean branch lengths rooted at the midpoint:"
  print $ toNewick $ measurableToPhyloTree tr
  let fnMeanTree = bnAnalysis ++ ".meantree"
  putStrLn $ "Save the mean tree to " <> fnMeanTree <> "."
  BL.writeFile fnMeanTree (toNewick $ measurableToPhyloTree tr)

  putStrLn "Root the trees at the midpoint of the mean tree."
  let outgroups = fst $ fromBipartition $ either error id $ bipartition tr
      trsRooted = map (first fromLength . either error id . outgroup outgroups "root" . first Length) trs
  putStrLn "Get the posterior means and the posterior covariance matrix."
  let pmR = getPosteriorMatrixRooted trsRooted
      (mu, sigma) = L.meanCov pmR
      (sigmaInv, (logSigmaDet, _)) = L.invlndet $ L.unSym sigma

  let fnData = bnAnalysis ++ ".data"
  putStrLn $ "Save the posterior means and covariances to " <> fnData <> "."
  encodeFile fnData (mu, L.toRows sigmaInv, logSigmaDet)
