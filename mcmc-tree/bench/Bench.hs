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

-- type Lens s t a b = forall f. Functor f => (a -> f b) -> s -> f t
-- type Lens' s a = forall f. Functor f => (a -> f a) -> s -> f s
-- lens :: (s -> a) -> (s -> b -> t) -> Lens s t a b
-- lens sa sbt afb s = sbt s <$> afb (sa s)

splitAt' :: Int -> [a] -> ([a], a, [a])
splitAt' i xs = (ls, head rs, tail rs)
  where
    (ls, rs) = splitAt i xs

assemble :: e -> a -> [Tree e a] -> [Tree e a] -> Tree e a -> Tree e a
assemble br lb ls rs c = Node br lb $ ls ++ (c : rs)

subTreeAt' :: Path -> Lens' (Tree e a) (Tree e a)
-- subTreeAt' p afb s = sbt s <$> afb (sa s)
--   where sbt t t' = let pos = goPathUnsafe p $ fromTree t in toTree $ pos {current = t'}
--         sa = getSubTreeUnsafe p
subTreeAt' pth f s = go s pth
  where
    go t [] = f t
    go (Node lb br ts) (x : xs) =
      let (ls, c, rs) = splitAt' x ts
       in assemble lb br ls rs <$> go c xs

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
    [ bgroup
        "lens"
        [ bench "change leaf Gn, lens with zipper" $ nf (changeLeaf (subTreeAt pthGn . root) "") tr,
          bench "change leaf Gn, optimized lens" $ nf (changeLeaf (subTreeAt' pthGn . root) "") tr,
          bench "change leaf Br, lens with zipper" $ nf (changeLeaf (subTreeAt pthBr . root) "") tr,
          bench "change leaf Br, optimized lens" $ nf (changeLeaf (subTreeAt' pthBr . root) "") tr
        ]
    ]
