{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Main
-- Description :  Simple tests on trees
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
module Main
  ( main,
  )
where

import Control.Lens
import Control.Monad
import ELynx.Tree
import Mcmc
import Mcmc.Tree
import Numeric.Log
import System.Random.MWC hiding (uniform)

type I = Tree Length (HeightLabel Name)

-- Birth death prior.
pr :: I -> Log Double
pr = birthDeath WithoutStem 1.0 1.0 1.0

-- No data.
lh :: I -> Log Double
lh _ = 1.0

-- Proposals on the tree.
proposals :: Show a => Tree e a -> Cycle I
proposals t =
  fromList $
    -- -- Pulley on the root node.
    -- pulleyUltrametric t 0.1 (PName "Tree root") (PWeight 5) Tune :
    -- Scale branches excluding the stem.
    [ subTreeAtE pth
        @~ slideNodeUltrametric 0.1 (PName $ "Tree node " ++ show lb) (PWeight 1) Tune
      | (pth, lb) <- itoList $ identify t,
        not (null pth),
        let s = t ^. subTreeAtE pth,
        not $ null $ forest s
    ]
      -- Scale trees of inner nodes excluding the root and the leaves.
      ++ [ subTreeAtE pth
             @~ scaleSubTreeUltrametric s 100 (PName $ "Tree node " ++ show lb) (PWeight 1) Tune
           | (pth, lb) <- itoList $ identify t,
             let s = t ^. subTreeAtE pth,
             not $ null pth,
             not $ null $ forest s
         ]

-- Get the height of the node at path. Useful to have a look at calibrated nodes.
getTreeNodeHeight :: Path -> I -> Double
getTreeNodeHeight p t = fromLength $ rootHeight $ t ^. subTreeAtE p

-- Monitor the height of all nodes.
monPs :: Tree e a -> [MonitorParameter I]
monPs t =
  [ getTreeNodeHeight pth >$< monitorDouble ("Node " ++ show lb)
    | (pth, lb) <- itoList $ identify t,
      let s = t ^. subTreeAtE pth,
      -- Path does not lead to a leaf.
      not $ null $ forest s
  ]

-- Monitor to standard output.
monStd :: Tree e a -> MonitorStdOut I
monStd t = monitorStdOut (monPs t) 1

monFile :: Tree e a -> MonitorFile I
monFile t = monitorFile "" (monPs t) 1

monTree :: MonitorFile I
monTree = monitorFile "-tree" [fromHeightTree >$< monitorTree "Tree"] 1

-- Combine the monitors.
mon :: Tree e a -> Monitor I
mon t = Monitor (monStd t) [monFile t, monTree] []

-- Number of burn in iterations.
nBurnIn :: Maybe Int
nBurnIn = Just 2000

-- Auto tuning period.
nAutoTune :: Maybe Int
nAutoTune = Just 100

-- Number of Metropolis-Hastings iterations after burn in.
nIter :: Int
nIter = 6000

main :: IO ()
main = do
  let t =
        toHeightTree $
          either error id $
            phyloToLengthTree $
              parseNewick Standard "(((a:1.0,b:1.0):1.0,c:2.0):1.0,d:3.0):0.0;"
  g <- create
  let s = force $ status "test" pr lh (proposals t) (mon t) t nBurnIn nAutoTune nIter g
  void $ mh s
