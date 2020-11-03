{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- |
-- Module      :  Bench
-- Description :  Benchmark common functions
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Tue Nov  3 14:47:08 2020.
module Main
  ( main,
  )
where

import Control.Lens
import Criterion.Main
import Data.Maybe
import ELynx.Tree
import Mcmc.Tree

--   let bf1 =
--         toTree . insertLabel "Bla"
--           . fromMaybe (error "Path does not lead to a leaf.")
--           . goPath pth
--           . fromTree
--   putStrLn $ "Change a leaf: " <> show (bf1 tr) <> "."
--   putStrLn "Benchmark change a leaf."
--   benchmark $ nf bf1 tr
--   let bf2 =
--         label . current
--           . fromMaybe (error "Path does not lead to a leaf.")
--           . goPath pth
--           . fromTree
--   putStrLn $ "Leaf to get: " <> show (bf2 tr) <> "."
--   putStrLn "Benchmark get a leaf."
--   benchmark $ nf bf2 tr
--   putStrLn "Benchmark calculation of prior."
--   let i = initWith tr
--       pr' = priorDistribution (getCalibrations tr) (getConstraints tr)
--   benchmark $ nf pr' i
--   (Just (mu, sigmaInvRows, logSigmaDet)) <- decodeFileStrict' fnData
--   let sigmaInv = L.fromRows sigmaInvRows
--       lh' = likelihoodFunction mu sigmaInv logSigmaDet
--   putStrLn "Benchmark calculation of likelihood."
--   benchmark $ nf lh' i

--   putStrLn "Benchmark identify."
--   benchmark $ nf identify tr

changeLeaf :: Lens' (Tree e a) a -> a -> Tree e a -> Tree e a
changeLeaf l x t = t & l .~ x

getPath :: (Show a, Eq a) => a -> Tree e a -> Path
getPath x t = fst $ fromMaybe err $ ifind (\_ n -> n == x) t
  where err = error $ show x <> " not found."

main :: IO ()
main = do
  tr <- oneTree "data/bench.tree"
  let nameGn = "Gnetum_montanum"
      pthGn = getPath nameGn tr
  let nameBr = "Brachypodium_distachyon"
      pthBr = getPath nameBr tr
  -- -- Some debugging.
  -- putStrLn $ "The path to \"Gnetum_montanum\" is: " <> show pth <> "."
  -- print $ toNewick $ measurableToPhyloTree tr
  -- let tr' = changeLeaf (subTreeAt' pth . root) "" tr
  -- print $ toNewick $ measurableToPhyloTree tr'
  defaultMain
    -- Optimized lenses are around 10 percent faster for this tree.
    [ bgroup
        "lens"
        [ bench "change leaf Gn, optimized lens" $ nf (changeLeaf (subTreeAt' pthGn . root) "") tr,
          bench "change leaf Br, optimized lens" $ nf (changeLeaf (subTreeAt' pthBr . root) "") tr
        ]
    ]
