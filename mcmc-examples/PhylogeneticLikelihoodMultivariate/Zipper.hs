-- |
-- Module      :  Zipper
-- Description :  Zippers on trees
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed Jul 22 14:17:59 2020.
module Zipper
  (
    TreePos (..),
    fromTree,
    toTree,
    parent,
    root,
    left,
    right,
    child,
    focus,
  )
where

import Data.Foldable
import ELynx.Data.Tree

type Content e a = (e, a)

fromContent :: Content e a -> Forest e a -> Tree e a
fromContent (br, lb) = Node br lb

toContent :: Tree e a -> Content e a
toContent (Node br lb _) = (br, lb)

data TreePos e a = Pos
  { -- | The currently selected tree.
    current :: Tree e a,
    -- | Forest to the left in reversed order.
    before :: Forest e a,
    -- | Forest to the right
    after :: Forest e a,
    -- | Finger to the selected tree
    parents :: [([Tree e a], Content e a, [Tree e a])]
  }
  deriving (Show, Eq)

fromTree :: Tree e a -> TreePos e a
fromTree t = Pos {current = t, before = [], after = [], parents = []}

toTree :: TreePos e a -> Tree e a
toTree = current . root

forestAt :: TreePos e a -> Forest e a
forestAt pos = foldl (flip (:)) (current pos : after pos) (before pos)

parent :: TreePos e a -> Maybe (TreePos e a)
parent pos = case parents pos of
  (ls, c, rs) : ps ->
    Just
      Pos
        { current = fromContent c $ forestAt pos,
          before = ls,
          after = rs,
          parents = ps
        }
  [] -> Nothing

root :: TreePos e a -> TreePos e a
root pos = maybe pos root (parent pos)

left :: TreePos e a -> Maybe (TreePos e a)
left pos =
  case before pos of
    t : ts ->
      Just
        pos
          { current = t,
            before = ts,
            after = current pos : after pos
          }
    [] -> Nothing

right :: TreePos e a -> Maybe (TreePos e a)
right pos =
  case after pos of
    t : ts ->
      Just
        pos
          { current = t,
            before = current pos : before pos,
            after = ts
          }
    [] -> Nothing

child :: TreePos e a -> Int -> Maybe (TreePos e a)
child pos n = case current pos of
  (Node _ _ ts)
    | null ts -> Nothing
    | length ts <= n -> Nothing
    | otherwise ->
      Just $
        Pos
          { current = head rs',
            before = reverse ls',
            after = tail rs',
            parents = (before pos, toContent (current pos), after pos) : parents pos
          }
    where
      (ls', rs') = splitAt n ts

focus :: TreePos e a -> [Int] -> Maybe (TreePos e a)
focus = foldlM child
