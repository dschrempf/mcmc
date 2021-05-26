{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

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
import Data.Aeson
import ELynx.Tree
import GHC.Generics
import Mcmc
import Mcmc.Tree
import System.Random.MWC hiding (uniform)

data I = I
  { _timeTree :: HeightTree Name,
    _rateTree :: Tree Length Name
  }
  deriving (Generic)

makeLenses ''I

instance ToJSON I

instance FromJSON I

-- Birth death prior, see Yang (2006), Figure 7.12. For lambda=2, mu=2, rho=0.1,
-- we expect a relatively linear distribution of inner node ages (a little bit
-- biased to younger ages).
pr :: PriorFunction I
pr (I t r) =
  product'
    [ birthDeath ConditionOnTimeOfMrca 2.0 2.0 0.1 $ fromHeightTree t,
      branchesWith withoutStem (exponential 1.0 . fromLength) r
    ]

-- The branches leading to the root are split. This Jacobian is necessary to
-- have unbiased proposals on the branches leading to the root.
jacobianRootBranch :: JacobianFunction I
jacobianRootBranch (I tTr rTr) = Exp $ log $ recip $ t1 * r1 + t2 * r2
  where
    (t1, t2) = case tTr of
      Node _ _ [l, r] -> (fromHeight $ nodeHeight $ label l, fromHeight $ nodeHeight $ label r)
      _ -> error "jacobianRootBranch: Time tree is not bifurcating."
    (r1, r2) = case rTr of
      Node _ _ [l, r] -> (fromLength $ branch l, fromLength $ branch r)
      _ -> error "jacobianRootBranch: Rate tree is not bifurcating."

-- Proposals on the ultrametric tree.
psT :: Tree e a -> [Proposal I]
psT t = map (liftProposalWith jacobianRootBranch timeTree) psWithJacobian
        ++ map (liftProposal timeTree) psWithoutJacobian
  where
    n = PName "Ultrametric tree"
    w = PWeight 1
    ps hd = slideNodesUltrametric t hd 0.1 n w Tune ++ scaleSubTreesUltrametric t hd 0.1 n w Tune
    psWithJacobian = pulleyUltrametric t 0.1 n (PWeight 5) Tune : ps (== 1)
    psWithoutJacobian = ps (>1)

-- Proposals on the unconstrained tree.
psR :: Tree e a -> [Proposal I]
psR t = map (liftProposalWith jacobianRootBranch rateTree) psWithJacobian
        ++ map (liftProposal rateTree) psWithoutJacobian
  where
    n = PName "Unconstrained tree"
    w = PWeight 1
    ps hd = scaleBranches t hd 0.1 n (PWeight 1) Tune ++ scaleSubTrees t hd 100 n (PWeight 1) Tune
    psWithJacobian = pulley 0.1 n (PWeight 5) Tune : ps (== 1)
    psWithoutJacobian = ps (> 1)

cc :: Tree e a -> Cycle I
cc t = cycleFromList $ psT t ++ psR t

-- Get the height of the node at path. Useful to have a look at calibrated nodes.
getNodeHeightT :: Path -> I -> Double
getNodeHeightT p x = fromHeight $ nodeHeight $ label $ x ^. timeTree . subTreeAtUnsafeL p

getBranchR :: Path -> I -> Double
getBranchR p x = fromLength $ branch $ x ^. rateTree . subTreeAtUnsafeL p

-- Monitor the height of all nodes.
monPs :: Tree e a -> [MonitorParameter I]
monPs t =
  [ getNodeHeightT pth >$< monitorDouble ("T Node" ++ show lb)
    | (pth, lb) <- itoList $ identify t,
      let s = t ^. subTreeAtUnsafeL pth,
      -- Path does not lead to a leaf.
      not $ null $ forest s
  ]
    ++ [ getBranchR pth >$< monitorDouble ("R Branch" ++ show lb)
         | (pth, lb) <- itoList $ identify t,
           not $ null pth
       ]

-- Monitor to standard output.
monStd :: Tree e a -> MonitorStdOut I
monStd t = monitorStdOut (monPs t) 10

monFile :: Tree e a -> MonitorFile I
monFile t = monitorFile "" (monPs t) 5

monTreeT :: MonitorFile I
monTreeT = monitorFile "ultrametric" [fromHeightTree . _timeTree >$< monitorTree "Tree"] 5

monTreeR :: MonitorFile I
monTreeR = monitorFile "unconstrained" [_rateTree >$< monitorTree "Tree"] 5

-- Combine the monitors.
mon :: Tree e a -> Monitor I
mon t = Monitor (monStd t) [monFile t, monTreeT, monTreeR] []

main :: IO ()
main = do
  let r =
        normalizeHeight $
          either error id $
            phyloToLengthTree $
              either error id $
                parseOneNewick Standard "(((a:1.0,b:1.0):1.0,c:2.0):1.0,(d:2.0,e:2.0):1.0):0.0;"
      t = either error id $ toHeightTreeUltrametric r
      cc' = cc t
  print r
  print t
  g <- create
  let mcmcS =
        Settings
          (AnalysisName "test-tree")
          (BurnInWithCustomAutoTuning [100, 110 .. 500])
          (Iterations 20000)
          Overwrite
          Sequential
          NoSave
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg pr noLikelihood cc' (mon t) TraceAuto (I t r) g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings 3 2
  -- a <- mc3 mc3S pr noLikelihood cc' (mon t) (I t r) g
  void $ mcmc mcmcS a
