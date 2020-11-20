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
import System.Random.MWC hiding (uniform)

type I = Tree Length (HeightLabel Name)

-- Birth death prior, see Yang (2006), Figure 7.12. For lambda=2, mu=2, rho=0.1,
-- we expect a relatively linear distribution of inner node ages (a little bit
-- biased to younger ages).
pr :: PriorFunction I
pr = birthDeath WithoutStem 2.0 2.0 0.1

-- Proposals on the tree.
cc :: Show a => Tree e a -> Cycle I
cc t =
  cycleFromList $
    -- -- Pulley on the root node.
    pulleyUltrametric t 0.1 (PName "Tree root") (PWeight 5) Tune :
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

burnIn :: BurnIn
burnIn = BurnInWithAutoTuning 2000 1000

iterations :: Int
iterations = 6000

main :: IO ()
main = do
  let t =
        toHeightTree $
          normalizeHeight $
            either error id $
              phyloToLengthTree $
                parseNewick Standard "(((a:1.0,b:1.0):1.0,c:2.0):1.0,(d:2.0,e:2.0):1.0):0.0;"
      cc' = cc t
  g <- create
  let s = Settings "test" burnIn iterations Overwrite NoSave Info
      a = mhg pr noLikelihood cc' (mon t) t g
  void $ mcmc s a
