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
  -- XXX
  -- product'
  --   [ birthDeath ConditionOnTimeOfMrca 2.0 2.0 0.1 $ fromHeightTree t,
  --     branchesWith withoutStem (exponential 0.1 . fromLength) r
  --   ]
  -- product'
  --   [branchesWith withoutStem (normal 1.0 1.0 . fromLength) r]
  1.0

-- Let the root branch measured in expected number of substitutions be
-- exponentially distributed.
lh :: LikelihoodFunction I
lh x = product $ dRoot (rootBranch x) : zipWith (\t r -> dOthers (fromLength t * fromLength r)) ts rs
  where
    getBranches (Node _ _ [l, r]) = tail (branches l) ++ tail (branches r)
    getBranches _ = error "getBranches: Root node is not bifurcating."
    ts = getBranches $ fromHeightTree $ _timeTree x
    rs = getBranches $ _rateTree x
    -- dRoot = exponential 1
    -- dOthers = exponential 1.0
    dRoot = gamma 1 (recip 1)
    dOthers = gamma 1 (recip 3)

-- The branches leading to the root are split. This Jacobian is necessary to
-- have unbiased proposals on the branches leading to the root.
rootBranch :: I -> Double
rootBranch (I tTr rTr) = t1 * r1 + t2 * r2
  where
    (t1, t2) = case fromHeightTree tTr of
      Node _ _ [l, r] -> (fromLength $ branch l, fromLength $ branch r)
      _ -> error "jacobianRootBranch: Time tree is not bifurcating."
    (r1, r2) = case rTr of
      Node _ _ [l, r] -> (fromLength $ branch l, fromLength $ branch r)
      _ -> error "jacobianRootBranch: Rate tree is not bifurcating."

jacobianRootBranch :: JacobianFunction I
jacobianRootBranch = Exp . log . recip . rootBranch

-- Proposals on the ultrametric tree.
psT :: Tree e a -> [Proposal I]
psT t =
  map (liftProposalWith jacobianRootBranch timeTree) psWithJacobian
    ++ map (liftProposal timeTree) psWithoutJacobian
  where
    w = PWeight 1
    ps hd n = slideNodesUltrametric t hd 0.5 n w Tune ++ scaleSubTreesUltrametric t hd 0.5 n w Tune
    nJ = PName "Time tree (Jac)"
    psWithJacobian = pulleyUltrametric t 0.5 nJ (PWeight 5) Tune : ps (== 1) nJ
    n0 = PName "Time tree"
    psWithoutJacobian = ps (> 1) n0

-- Proposals on the unconstrained tree.
psR :: Tree e a -> [Proposal I]
psR t =
  map (liftProposalWith jacobianRootBranch rateTree) psWithJacobian
    ++ map (liftProposal rateTree) psWithoutJacobian
  where
    w = PWeight 1
    ps hd n = scaleBranches t hd 0.1 n w Tune ++ scaleSubTrees t hd 100 n w Tune
    nJ = PName "Rate tree (Jac)"
    psWithJacobian = pulley 0.1 nJ (PWeight 5) Tune : ps (== 1) nJ
    n0 = PName "Rate tree"
    psWithoutJacobian = ps (> 1) n0

cc :: Tree e a -> Cycle I
-- XXX.
-- cc t = cycleFromList $ psT t ++ psR t
cc t = cycleFromList $ psR t

-- Get the height of the node at path. Useful to have a look at calibrated nodes.
getTimeTreeNodeHeight :: Path -> I -> Double
getTimeTreeNodeHeight p x = fromHeight $ nodeHeight $ label $ x ^. timeTree . subTreeAtUnsafeL p

getTimeTreeBranch :: Path -> I -> Double
getTimeTreeBranch p x = fromLength $ branch $ fromHeightTree (x ^. timeTree) ^. subTreeAtUnsafeL p

getRateTreeBranch :: Path -> I -> Double
getRateTreeBranch p x = fromLength $ branch $ x ^. rateTree . subTreeAtUnsafeL p

-- Monitor the height of all nodes.
monPs :: Tree e a -> [MonitorParameter I]
monPs t =
  (rootBranch >$< monitorDouble "Root Branch") :
  [ (\x -> getTimeTreeBranch pth x * getRateTreeBranch pth x) >$< monitorDouble ("T*R Branch" ++ show lb)
    | (pth, lb) <- itoList $ identify t,
      length pth > 1
  ]
    ++ [ getTimeTreeNodeHeight pth >$< monitorDouble ("T Node" ++ show lb)
         | (pth, lb) <- itoList $ identify t,
           let s = t ^. subTreeAtUnsafeL pth,
           -- Path does not lead to a leaf.
           not $ null $ forest s
       ]
    ++ [ getTimeTreeBranch pth >$< monitorDouble ("T Branch" ++ show lb)
         | (pth, lb) <- itoList $ identify t,
           not $ null pth
       ]
    ++ [ getRateTreeBranch pth >$< monitorDouble ("R Branch" ++ show lb)
         | (pth, lb) <- itoList $ identify t,
           not $ null pth
       ]

-- Monitor to standard output.
monStd :: Tree e a -> MonitorStdOut I
monStd t = monitorStdOut (take 4 $ monPs t) 100

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
        either error id $
          phyloToLengthTree $
            either error id $
              parseOneNewick Standard "(((a:1.0,b:1.0):1.0,(c:1.0,d:1.0):1.0):1.0,(e:1.0,f:1.0):2.0):0.0;"
      t = either error id $ toHeightTreeUltrametric $ normalizeHeight r
  print r
  print t
  g <- create
  let mcmcS =
        Settings
          (AnalysisName "test-tree")
          (BurnInWithCustomAutoTuning [100, 110 .. 500])
          (Iterations 30000)
          Overwrite
          Sequential
          NoSave
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg pr lh (cc t) (mon t) TraceAuto (I t r) g
  -- -- Metropolic-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings 3 2
  -- a <- mc3 mc3S pr noLikelihood cc' (mon t) (I t r) g
  void $ mcmc mcmcS a
