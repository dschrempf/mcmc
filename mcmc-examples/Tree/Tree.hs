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
--
-- Perform an MCMC run on a pair of trees. This example solves the phylogenetic
-- problem of dating (estimating node ages) a tree.
module Main
  ( main,
  )
where

import Control.Lens
import Control.Monad
import Data.Aeson
import Data.Bifunctor
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

-- -- Birth death prior, see Yang (2006), Figure 7.12. For lambda=2, mu=2, rho=0.1,
-- -- we expect a relatively linear distribution of inner node ages (a little bit
-- -- biased to younger ages). A hyper prior for the birth and death rates would be
-- -- required.
-- pr :: PriorFunction I
-- pr (I t r) =
--   product'
--     [ birthDeath ConditionOnTimeOfMrca 2.0 2.0 0.1 $ fromHeightTree t,
--       branchesWith withoutStem (normal 1.0 0.4 . fromLength) r
--     ]

pr :: PriorFunction I
pr (I _ r) =
  product'
    [branchesWith withoutStem (normal 1.0 0.3 . fromLength) r]

-- The branches measured in expected number of substitutions are determining the likelihood.
lh :: LikelihoodFunction I
lh x = product $ dRoot (rootBranch x) : zipWith (\t r -> dOthers (fromLength t * fromLength r)) ts rs
  where
    getBranches (Node _ _ [l, r]) = tail (branches l) ++ tail (branches r)
    getBranches _ = error "getBranches: Root node is not bifurcating."
    ts = getBranches $ fromHeightTree $ _timeTree x
    rs = getBranches $ _rateTree x
    dRoot = normal 1 0.1
    dOthers = normal (1 / 3) 0.1

-- The root splits the branch of the unrooted tree into two branches. This
-- function retrieves the root branch measured in expected number of
-- substitutions.
rootBranch :: I -> Double
rootBranch (I tTr rTr) = t1 * r1 + t2 * r2
  where
    (t1, t2) = case fromHeightTree tTr of
      Node _ _ [l, r] -> (fromLength $ branch l, fromLength $ branch r)
      _ -> error "rootBranch: Time tree is not bifurcating."
    (r1, r2) = case rTr of
      Node _ _ [l, r] -> (fromLength $ branch l, fromLength $ branch r)
      _ -> error "rootBranch: Rate tree is not bifurcating."

-- This Jacobian is necessary to have unbiased proposals on the branches leading
-- to the root.
jacobianRootBranch :: JacobianFunction I
jacobianRootBranch = Exp . log . recip . rootBranch

-- Proposals on the time tree.
psT :: Tree e a -> [Proposal I]
psT t =
  map (liftProposalWith jacobianRootBranch timeTree) psWithJacobian
    ++ map (liftProposal timeTree) psWithoutJacobian
  where
    w = PWeight 1
    ps hd n = slideNodesUltrametric t hd 0.5 n w Tune ++ scaleSubTreesUltrametric t hd 0.5 n w Tune
    nR = PName "Time tree (root branches)"
    psWithJacobian = pulleyUltrametric t 0.5 nR w Tune : ps (== 1) nR
    nO = PName "Time tree"
    psWithoutJacobian = ps (> 1) nO

-- Proposals on the rate tree.
psR :: Tree e a -> [Proposal I]
psR t =
  map (liftProposalWith jacobianRootBranch rateTree) psAtRoot
    ++ map (liftProposal rateTree) psOthers
  where
    w = PWeight 1
    ps hd n = scaleBranches t hd 5.0 n w Tune ++ scaleSubTrees t hd 100 n w Tune
    nR = PName "Rate tree (root branches)"
    psAtRoot = pulley 0.5 nR w Tune : ps (==1) nR
    nO = PName "Rate tree"
    psOthers = scaleBranches t (> 1) 5.0 nO w Tune ++ scaleSubTrees t (> 1) 100 nO w Tune

-- A contrary proposal on the time and rate trees.
psContra :: Tree e a -> [Proposal I]
psContra t =
  map (liftProposalWith jacobianRootBranch timeRateTreesL) psAtRoot
    ++ map (liftProposal timeRateTreesL) psOthers
  where
    psAtRoot = scaleSubTreesContrarily t (== 1) 0.01 (PName "Trees contra") (PWeight 3) Tune
    psOthers = scaleSubTreesContrarily t (> 1) 0.01 (PName "Trees contra") (PWeight 3) Tune
    -- Lens for a contrary proposal on the trees.
    timeRateTreesL :: Lens' I (HeightTree Name, Tree Length Name)
    timeRateTreesL =
      lens
        (\x -> (x ^. timeTree, x ^. rateTree))
        (\x (tTr, rTr) -> x {_timeTree = tTr, _rateTree = rTr})

-- The cycle includes proposals on both trees.
cc :: Tree e a -> Cycle I
cc t = cycleFromList $ psT t ++ psR t ++ psContra t

getTimeTreeNodeHeight :: Path -> I -> Double
getTimeTreeNodeHeight p x = fromHeight $ nodeHeight $ label $ x ^. timeTree . subTreeAtUnsafeL p

getTimeTreeBranch :: Path -> I -> Double
getTimeTreeBranch p x = fromLength $ branch $ fromHeightTree (x ^. timeTree) ^. subTreeAtUnsafeL p

getRateTreeBranch :: Path -> I -> Double
getRateTreeBranch p x = fromLength $ branch $ x ^. rateTree . subTreeAtUnsafeL p

-- Useful monitors.
monPs :: Tree e a -> [MonitorParameter I]
monPs t =
  -- Branches measured in expected number of substitutions.
  (rootBranch >$< monitorDouble "Root Branch") :
  [ (\x -> getTimeTreeBranch pth x * getRateTreeBranch pth x) >$< monitorDouble ("T*R Branch" ++ show lb)
    | (pth, lb) <- itoList $ identify t,
      length pth > 1
  ]
    -- Node heights of the time tree.
    ++ [ getTimeTreeNodeHeight pth >$< monitorDouble ("T Node" ++ show lb)
         | (pth, lb) <- itoList $ identify t,
           let s = t ^. subTreeAtUnsafeL pth,
           -- Path does not lead to a leaf.
           not $ null $ forest s
       ]
    -- Branch lengths of the time tree.
    ++ [ getTimeTreeBranch pth >$< monitorDouble ("T Branch" ++ show lb)
         | (pth, lb) <- itoList $ identify t,
           not $ null pth
       ]
    -- Branch lengths of the rate tree.
    ++ [ getRateTreeBranch pth >$< monitorDouble ("R Branch" ++ show lb)
         | (pth, lb) <- itoList $ identify t,
           not $ null pth
       ]

monStd :: Tree e a -> MonitorStdOut I
monStd t = monitorStdOut (take 4 $ monPs t) 100

monFile :: Tree e a -> MonitorFile I
monFile t = monitorFile "" (monPs t) 5

monTreeT :: MonitorFile I
monTreeT = monitorFile "ultrametric" [fromHeightTree . _timeTree >$< monitorTree "Tree"] 5

monTreeR :: MonitorFile I
monTreeR = monitorFile "unconstrained" [_rateTree >$< monitorTree "Tree"] 5

mon :: Tree e a -> Monitor I
mon t = Monitor (monStd t) [monFile t, monTreeT, monTreeR] []

main :: IO ()
main = do
  let -- The tree to be dated.
      tree = "(((a:1.0,b:1.0):1.0,(c:1.0,d:1.0):1.0):1.0,(e:1.0,f:1.0):2.0):0.0;"
      r' =
        either error id $
          phyloToLengthTree $
            either error id $
              parseOneNewick Standard tree
      t = either error id $ toHeightTreeUltrametric $ normalizeHeight r'
      -- Set the initial rates to 1.0.
      r = first (const 1.0) r'
  print r
  print t
  g <- create
  let mcmcS =
        Settings
          (AnalysisName "test-tree")
          (BurnInWithCustomAutoTuning [100, 110 .. 500])
          (Iterations 40000)
          Overwrite
          Sequential
          NoSave
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg pr lh (cc t) (mon t) TraceAuto (I t r) g
  -- -- Metropolis-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings (NChains 3) (SwapPeriod 2) (NSwaps 1)
  -- a <- mc3 mc3S pr lh (cc t) (mon t) TraceAuto (I t r) g
  void $ mcmc mcmcS a
