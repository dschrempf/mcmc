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
-- Perform an MCMC run to date a phylogenetic tree (estimate the node ages of a
-- phylogenetic tree). This example involves a pair of trees. The first tree is
-- the "time tree" with branch lengths in units proportional to real time. The
-- second tree is the "rate tree", which has the same topology as the time tree.
-- The branches of the rate tree describe the rate of evolution at the very same
-- branch of the time tree.
--
-- For a detailed description, please see
-- https://github.com/dschrempf/mcmc-dating/blob/master/src/Definitions.hs.
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

-- The state space.
data I = I
  { _timeBirthRate :: Double,
    _timeDeathRate :: Double,
    _timeTree :: HeightTree Name,
    _rateMean :: Double,
    _rateVariance :: Double,
    _rateTree :: Tree Length Name
  }
  deriving (Generic)

makeLenses ''I

instance ToJSON I

instance FromJSON I

-- The prior functions.
pr :: PriorFunction I
pr (I l m t mu va r) =
  product'
    [ exponential 1 l,
      exponential 1 m,
      birthDeath ConditionOnTimeOfMrca l m 1.0 $ fromHeightTree t,
      exponential 1 mu,
      exponential 1 va,
      uncorrelatedGamma withoutStem 1 va r
    ]

-- The branches measured in expected number of substitutions (substitutions =
-- time * rate) are determining the likelihood.
lh :: LikelihoodFunction I
lh x@(I _ _ t mu _ r) =
  product' $
    dRoot (rootBranch x) : zipWith (\tBr rBr -> dOthers (fromLength tBr * mu * fromLength rBr)) ts rs
  where
    getBranches (Node _ _ [lTr, rTr]) = tail (branches lTr) ++ tail (branches rTr)
    getBranches _ = error "getBranches: Root node is not bifurcating."
    ts = getBranches $ fromHeightTree t
    rs = getBranches r
    dRoot = normal 1.0 0.1
    dOthers = normal (1 / 3) 0.1

-- The root splits the branch of the unrooted tree into two branches. This
-- function retrieves the root branch measured in expected number of
-- substitutions.
rootBranch :: I -> Double
rootBranch (I _ _ tTr mu _ rTr) = mu * (t1 * r1 + t2 * r2)
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
  map (liftProposalWith jacobianRootBranch timeTree) psAtRoot
    ++ map (liftProposal timeTree) psOthers
  where
    nR = PName "Time tree [R]"
    nO = PName "Time tree [O]"
    ps hd n =
      slideNodesUltrametric t hd 0.5 n (PWeight 3) Tune
        ++ scaleSubTreesUltrametric t hd 0.5 n 1 Tune
    psAtRoot = pulleyUltrametric t 0.5 nR (PWeight 6) Tune : ps (== 1) nR
    psOthers = ps (> 1) nO

-- Proposals on the rate tree.
psR :: Tree e a -> [Proposal I]
psR t =
  map (liftProposalWith jacobianRootBranch rateTree) psAtRoot
    ++ map (liftProposal rateTree) psOthers
  where
    nR = PName "Rate tree [R]"
    nO = PName "Rate tree [O]"
    ps hd n =
      scaleBranches t hd 5.0 n (PWeight 3) Tune
        ++ scaleSubTrees t hd 100 n 1 Tune
    psAtRoot = pulley 0.5 nR (PWeight 6) Tune : ps (== 1) nR
    psOthers = ps (> 1) nO

-- A contrary proposal on the time and rate trees.
psContra :: Tree e a -> [Proposal I]
psContra t =
  map (liftProposalWith jacobianRootBranch timeRateTreesL) psAtRoot
    ++ map (liftProposal timeRateTreesL) psOthers
  where
    -- Lens for a contrary proposal on the trees.
    timeRateTreesL :: Lens' I (HeightTree Name, Tree Length Name)
    timeRateTreesL =
      lens
        (\x -> (x ^. timeTree, x ^. rateTree))
        (\x (tTr, rTr) -> x {_timeTree = tTr, _rateTree = rTr})
    w = PWeight 3
    nR = PName "Trees contra [R]"
    nO = PName "Trees contra [O]"
    ps hd n = scaleSubTreesContrarily t hd 0.01 n w Tune
    psAtRoot = ps (== 1) nR
    psOthers = ps (> 1) nO

-- The cycle includes proposals on both trees.
cc :: Tree e a -> Cycle I
cc t =
  cycleFromList $
    liftProposal timeBirthRate (scaleUnbiased 3.0 (PName "Birth rate") (PWeight 20) Tune) :
    liftProposal timeDeathRate (scaleUnbiased 3.0 (PName "Death rate") (PWeight 20) Tune) :
    liftProposal rateMean (scaleUnbiased 3.0 (PName "Rate mean") (PWeight 20) Tune) :
    liftProposal rateVariance (scaleUnbiased 3.0 (PName "Rate variance") (PWeight 20) Tune) :
    psT t ++ psR t ++ psContra t

getTimeTreeNodeHeight :: Path -> I -> Double
getTimeTreeNodeHeight p x = fromHeight $ nodeHeight $ label $ x ^. timeTree . subTreeAtUnsafeL p

getTimeTreeBranch :: Path -> I -> Double
getTimeTreeBranch p x = fromLength $ branch $ fromHeightTree (x ^. timeTree) ^. subTreeAtUnsafeL p

getRateTreeBranch :: Path -> I -> Double
getRateTreeBranch p x = fromLength $ branch $ x ^. rateTree . subTreeAtUnsafeL p

-- Useful monitors.
monPs :: Tree e a -> [MonitorParameter I]
monPs t =
  (view timeBirthRate >$< monitorDouble "Birth rate") :
  (view timeDeathRate >$< monitorDouble "Death rate") :
  (view rateMean >$< monitorDouble "Rate mean") :
  (view rateVariance >$< monitorDouble "Rate variance") :
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
      r = setStem 0 $ first (const 1.0) r'
  print r
  print t
  g <- create
  let mcmcS =
        Settings
          (AnalysisName "tree")
          (BurnInWithCustomAutoTuning [100, 110 .. 500])
          (Iterations 20000)
          Overwrite
          Sequential
          NoSave
          LogStdOutAndFile
          Info
  -- Metropolis-Hastings-Green algorithm.
  a <- mhg pr lh (cc t) (mon t) TraceAuto (I 1.0 1.0 t 1.0 1.0 r) g
  -- -- Metropolis-coupled Markov chain Monte Carlo algorithm.
  -- let mc3S = MC3Settings (NChains 3) (SwapPeriod 2) (NSwaps 1)
  -- a <- mc3 mc3S pr lh (cc t) (mon t) TraceAuto (I t r) g
  void $ mcmc mcmcS a
