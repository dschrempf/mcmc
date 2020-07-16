{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TupleSections #-}

-- |
-- Module     : AdjacencyIntMap
-- Copyright  : (c) Andrey Mokhov 2016-2019
-- License    : MIT (see the file LICENSE)
-- Maintainer : andrey.mokhov@gmail.com
-- Stability  : experimental
--
-- __Alga__ is a library for algebraic construction and manipulation of graphs
-- in Haskell. See <https://github.com/snowleopard/alga-paper this paper> for the
-- motivation behind the library, the underlying theory, and implementation details.
--
-- This module defines the 'AdjacencyMap' data type for edge-labelled graphs, as
-- well as associated operations and algorithms. 'AdjacencyMap' is an instance
-- of the 'C.Graph' type class, which can be used for polymorphic graph
-- construction and manipulation.
module AdjacencyIntMap
  ( -- * Data structure
    AdjacencyIntMap,
    adjacencyIntMap,

    -- * Basic graph construction primitives
    empty,
    vertex,
    edge,
    (-<),
    (>-),
    overlay,
    connect,
    vertices,
    edges,
    overlays,
    fromAdjacencyIntMaps,

    -- * Relations on graphs
    isSubgraphOf,

    -- * Graph properties
    isEmpty,
    hasVertex,
    hasEdge,
    edgeLabel,
    vertexCount,
    edgeCount,
    vertexList,
    edgeList,
    vertexSet,
    edgeSet,
    preSet,
    postSet,
    skeleton,

    -- * Graph transformation
    removeVertex,
    removeEdge,
    replaceVertex,
    replaceEdge,
    transpose,
    gmap,
    emap,
    induce,
    -- induceJust,

    -- * Relational operations
    closure,
    reflexiveClosure,
    symmetricClosure,
    transitiveClosure,

    -- * Miscellaneous
    consistent,
  )
where

import qualified Algebra.Graph.AdjacencyIntMap as AM
import Algebra.Graph.Label
import Control.DeepSeq
import Data.IntMap (IntMap)
import qualified Data.IntMap.Strict as M
import Data.Maybe
import Data.Monoid (Sum (..))
import Data.IntSet (IntSet, (\\))
import Data.Set (Set)
import qualified Data.Set as S
import qualified Data.IntSet as IS
import GHC.Generics

-- | Edge-labelled graphs, where the type variable @e@ stands for edge labels.
-- For example, 'AdjacencyMap' @Bool@ @a@ is isomorphic to unlabelled graphs
-- defined in the top-level module "Algebra.Graph.AdjacencyMap", where @False@
-- and @True@ denote the lack of and the existence of an unlabelled edge,
-- respectively.
newtype AdjacencyIntMap e
  = AM
      { -- | The /adjacency map/ of an edge-labelled graph: each vertex is
        -- associated with a map from its direct successors to the corresponding
        -- edge labels.
        adjacencyIntMap :: IntMap (IntMap e)
      }
  deriving (Eq, Generic, NFData)

instance (Ord e, Show e) => Show (AdjacencyIntMap e) where
  showsPrec p lam@(AM m)
    | IS.null vs = showString "empty"
    | null es = showParen (p > 10) $ vshow vs
    | vs == used = showParen (p > 10) $ eshow es
    | otherwise =
      showParen (p > 10) $
        showString "overlay (" . vshow (vs \\ used)
          . showString ") ("
          . eshow es
          . showString ")"
    where
      vs = vertexSet lam
      es = edgeList lam
      used = referredToVertexSet m
      vshow vs' = case IS.toAscList vs' of
        [x] -> showString "vertex " . showsPrec 11 x
        xs -> showString "vertices " . showsPrec 11 xs
      eshow es' = case es' of
        [(e, x, y)] ->
          showString "edge " . showsPrec 11 e
            . showString " "
            . showsPrec 11 x
            . showString " "
            . showsPrec 11 y
        xs -> showString "edges " . showsPrec 11 xs

instance (Ord e, Monoid e) => Ord (AdjacencyIntMap e) where
  compare x y =
    mconcat
      [ compare (vertexCount x) (vertexCount y),
        compare (vertexSet x) (vertexSet y),
        compare (edgeCount x) (edgeCount y),
        compare (eSet x) (eSet y),
        cmp
      ]
    where
      eSet = S.map (\(_, x', y') -> (x', y')) . edgeSet
      cmp
        | x == y = EQ
        | overlays [x, y] == y = LT
        | otherwise = compare x y

-- | __Note:__ this does not satisfy the usual ring laws; see 'AdjacencyMap'
-- for more details.
instance (Eq e, Dioid e) => Num (AdjacencyIntMap e) where
  fromInteger = vertex . fromInteger
  (+) = overlay
  (*) = connect mempty
  signum = const empty
  abs = id
  negate = id

-- | Construct the /empty graph/.
-- Complexity: /O(1)/ time and memory.
--
-- @
-- 'isEmpty'     empty == True
-- 'hasVertex' x empty == False
-- 'vertexCount' empty == 0
-- 'edgeCount'   empty == 0
-- @
empty :: AdjacencyIntMap e
empty = AM M.empty

-- | Construct the graph comprising /a single isolated vertex/.
-- Complexity: /O(1)/ time and memory.
--
-- @
-- 'isEmpty'     (vertex x) == False
-- 'hasVertex' x (vertex y) == (x == y)
-- 'vertexCount' (vertex x) == 1
-- 'edgeCount'   (vertex x) == 0
-- @
vertex :: Int -> AdjacencyIntMap e
vertex x = AM $ M.singleton x M.empty

-- | Construct the graph comprising /a single edge/.
-- Complexity: /O(1)/ time, memory.
--
-- @
-- edge e    x y              == 'connect' e ('vertex' x) ('vertex' y)
-- edge 'zero' x y              == 'vertices' [x,y]
-- 'hasEdge'   x y (edge e x y) == (e /= 'zero')
-- 'edgeLabel' x y (edge e x y) == e
-- 'edgeCount'     (edge e x y) == if e == 'zero' then 0 else 1
-- 'vertexCount'   (edge e 1 1) == 1
-- 'vertexCount'   (edge e 1 2) == 2
-- @
edge :: (Eq e, Monoid e) => e -> Int -> Int -> AdjacencyIntMap e
edge e x y
  | e == zero = vertices [x, y]
  | x == y = AM $ M.singleton x (M.singleton x e)
  | otherwise = AM $ M.fromList [(x, M.singleton y e), (y, M.empty)]

-- | The left-hand part of a convenient ternary-ish operator @x-\<e\>-y@ for
-- creating labelled edges.
--
-- @
-- x -\<e\>- y == 'edge' e x y
-- @
(-<) :: a -> e -> (a, e)
g -< e = (g, e)

-- | The right-hand part of a convenient ternary-ish operator @x-\<e\>-y@ for
-- creating labelled edges.
--
-- @
-- x -\<e\>- y == 'edge' e x y
-- @
(>-) :: (Eq e, Monoid e) => (Int, e) -> Int -> AdjacencyIntMap e
(x, e) >- y = edge e x y

infixl 5 -<

infixl 5 >-

-- | /Overlay/ two graphs. This is a commutative, associative and idempotent
-- operation with the identity 'empty'.
-- Complexity: /O((n + m) * log(n))/ time and /O(n + m)/ memory.
--
-- @
-- 'isEmpty'     (overlay x y) == 'isEmpty'   x   && 'isEmpty'   y
-- 'hasVertex' z (overlay x y) == 'hasVertex' z x || 'hasVertex' z y
-- 'vertexCount' (overlay x y) >= 'vertexCount' x
-- 'vertexCount' (overlay x y) <= 'vertexCount' x + 'vertexCount' y
-- 'edgeCount'   (overlay x y) >= 'edgeCount' x
-- 'edgeCount'   (overlay x y) <= 'edgeCount' x   + 'edgeCount' y
-- 'vertexCount' (overlay 1 2) == 2
-- 'edgeCount'   (overlay 1 2) == 0
-- @
--
-- Note: 'overlay' composes edges in parallel using the operator '<+>' with
-- 'zero' acting as the identity:
--
-- @
-- 'edgeLabel' x y $ overlay ('edge' e x y) ('edge' 'zero' x y) == e
-- 'edgeLabel' x y $ overlay ('edge' e x y) ('edge' f    x y) == e '<+>' f
-- @
--
-- Furthermore, when applied to transitive graphs, 'overlay' composes edges in
-- sequence using the operator '<.>' with 'one' acting as the identity:
--
-- @
-- 'edgeLabel' x z $ 'transitiveClosure' (overlay ('edge' e x y) ('edge' 'one' y z)) == e
-- 'edgeLabel' x z $ 'transitiveClosure' (overlay ('edge' e x y) ('edge' f   y z)) == e '<.>' f
-- @
overlay :: (Eq e, Monoid e) => AdjacencyIntMap e -> AdjacencyIntMap e -> AdjacencyIntMap e
overlay (AM x) (AM y) = AM $ M.unionWith nonZeroUnion x y

-- Union maps, removing zero elements from the result.
nonZeroUnion :: (Eq e, Monoid e) => IntMap e -> IntMap e -> IntMap e
nonZeroUnion x y = M.filter (/= zero) $ M.unionWith mappend x y

-- Drop all edges with zero labels.
trimZeroes :: (Eq e, Monoid e) => IntMap (IntMap e) -> IntMap (IntMap e)
trimZeroes = M.map (M.filter (/= zero))

-- | /Connect/ two graphs with edges labelled by a given label. When applied to
-- the same labels, this is an associative operation with the identity 'empty',
-- which distributes over 'overlay' and obeys the decomposition axiom.
-- Complexity: /O((n + m) * log(n))/ time and /O(n + m)/ memory. Note that the
-- number of edges in the resulting graph is quadratic with respect to the
-- number of vertices of the arguments: /m = O(m1 + m2 + n1 * n2)/.
--
-- @
-- 'isEmpty'     (connect e x y) == 'isEmpty'   x   && 'isEmpty'   y
-- 'hasVertex' z (connect e x y) == 'hasVertex' z x || 'hasVertex' z y
-- 'vertexCount' (connect e x y) >= 'vertexCount' x
-- 'vertexCount' (connect e x y) <= 'vertexCount' x + 'vertexCount' y
-- 'edgeCount'   (connect e x y) <= 'vertexCount' x * 'vertexCount' y + 'edgeCount' x + 'edgeCount' y
-- 'vertexCount' (connect e 1 2) == 2
-- 'edgeCount'   (connect e 1 2) == if e == 'zero' then 0 else 1
-- @
connect :: (Eq e, Monoid e) => e -> AdjacencyIntMap e -> AdjacencyIntMap e -> AdjacencyIntMap e
connect e (AM x) (AM y)
  | e == mempty = overlay (AM x) (AM y)
  | otherwise =
    AM $ M.unionsWith nonZeroUnion $
      x : y
        : [M.fromSet (const targets) (M.keysSet x)]
  where
    targets = M.fromSet (const e) (M.keysSet y)

-- | Construct the graph comprising a given list of isolated vertices.
-- Complexity: /O(L * log(L))/ time and /O(L)/ memory, where /L/ is the length
-- of the given list.
--
-- @
-- vertices []            == 'empty'
-- vertices [x]           == 'vertex' x
-- 'hasVertex' x . vertices == 'elem' x
-- 'vertexCount' . vertices == 'length' . 'Data.List.nub'
-- 'vertexSet'   . vertices == Set.'Set.fromList'
-- @
vertices :: [Int] -> AdjacencyIntMap e
vertices = AM . M.fromList . map (,M.empty)

-- | Construct the graph from a list of edges.
-- Complexity: /O((n + m) * log(n))/ time and /O(n + m)/ memory.
--
-- @
-- edges []        == 'empty'
-- edges [(e,x,y)] == 'edge' e x y
-- edges           == 'overlays' . 'map' (\\(e, x, y) -> 'edge' e x y)
-- @
edges :: (Eq e, Monoid e) => [(e, Int, Int)] -> AdjacencyIntMap e
edges es = fromAdjacencyIntMaps [(x, M.singleton y e) | (e, x, y) <- es]

-- | Overlay a given list of graphs.
-- Complexity: /O((n + m) * log(n))/ time and /O(n + m)/ memory.
--
-- @
-- overlays []        == 'empty'
-- overlays [x]       == x
-- overlays [x,y]     == 'overlay' x y
-- overlays           == 'foldr' 'overlay' 'empty'
-- 'isEmpty' . overlays == 'all' 'isEmpty'
-- @
overlays :: (Eq e, Monoid e) => [AdjacencyIntMap e] -> AdjacencyIntMap e
overlays = AM . M.unionsWith nonZeroUnion . map adjacencyIntMap

-- | Construct a graph from a list of adjacency sets.
-- Complexity: /O((n + m) * log(n))/ time and /O(n + m)/ memory.
--
-- @
-- fromAdjacencyIntMaps []                                  == 'empty'
-- fromAdjacencyIntMaps [(x, M.'M.empty')]                    == 'vertex' x
-- fromAdjacencyIntMaps [(x, M.'M.singleton' y e)]            == if e == 'zero' then 'vertices' [x,y] else 'edge' e x y
-- 'overlay' (fromAdjacencyIntMaps xs) (fromAdjacencyIntMaps ys) == fromAdjacencyIntMaps (xs '++' ys)
-- @
fromAdjacencyIntMaps :: (Eq e, Monoid e) => [(Int, IntMap e)] -> AdjacencyIntMap e
fromAdjacencyIntMaps xs = AM $ trimZeroes $ M.unionWith mappend vs es
  where
    vs = M.fromSet (const M.empty) . IS.unions $ map (M.keysSet . snd) xs
    es = M.fromListWith (M.unionWith mappend) xs

-- | The 'isSubgraphOf' function takes two graphs and returns 'True' if the
-- first graph is a /subgraph/ of the second.
-- Complexity: /O(s + m * log(m))/ time. Note that the number of edges /m/ of a
-- graph can be quadratic with respect to the expression size /s/.
--
-- @
-- isSubgraphOf 'empty'      x     ==  True
-- isSubgraphOf ('vertex' x) 'empty' ==  False
-- isSubgraphOf x y              ==> x <= y
-- @
isSubgraphOf :: (Eq e, Monoid e) => AdjacencyIntMap e -> AdjacencyIntMap e -> Bool
isSubgraphOf (AM x) (AM y) = M.isSubmapOfBy (M.isSubmapOfBy le) x y
  where
    le x' y' = mappend x' y' == y'

-- | Check if a graph is empty.
-- Complexity: /O(1)/ time.
--
-- @
-- isEmpty 'empty'                         == True
-- isEmpty ('overlay' 'empty' 'empty')         == True
-- isEmpty ('vertex' x)                    == False
-- isEmpty ('removeVertex' x $ 'vertex' x)   == True
-- isEmpty ('removeEdge' x y $ 'edge' e x y) == False
-- @
isEmpty :: AdjacencyIntMap e -> Bool
isEmpty = M.null . adjacencyIntMap

-- | Check if a graph contains a given vertex.
-- Complexity: /O(log(n))/ time.
--
-- @
-- hasVertex x 'empty'            == False
-- hasVertex x ('vertex' y)       == (x == y)
-- hasVertex x . 'removeVertex' x == 'const' False
-- @
hasVertex :: Int -> AdjacencyIntMap e -> Bool
hasVertex x = M.member x . adjacencyIntMap

-- | Check if a graph contains a given edge.
-- Complexity: /O(log(n))/ time.
--
-- @
-- hasEdge x y 'empty'            == False
-- hasEdge x y ('vertex' z)       == False
-- hasEdge x y ('edge' e x y)     == (e /= 'zero')
-- hasEdge x y . 'removeEdge' x y == 'const' False
-- hasEdge x y                  == 'not' . 'null' . 'filter' (\\(_,ex,ey) -> ex == x && ey == y) . 'edgeList'
-- @
hasEdge :: Int -> Int -> AdjacencyIntMap e -> Bool
hasEdge x y (AM m) = maybe False (M.member y) (M.lookup x m)

-- | Extract the label of a specified edge in a graph.
-- Complexity: /O(log(n))/ time.
--
-- @
-- edgeLabel x y 'empty'         == 'zero'
-- edgeLabel x y ('vertex' z)    == 'zero'
-- edgeLabel x y ('edge' e x y)  == e
-- edgeLabel s t ('overlay' x y) == edgeLabel s t x <+> edgeLabel s t y
-- @
edgeLabel :: Monoid e => Int -> Int -> AdjacencyIntMap e -> e
edgeLabel x y (AM m) = fromMaybe zero (M.lookup x m >>= M.lookup y)

-- | The number of vertices in a graph.
-- Complexity: /O(1)/ time.
--
-- @
-- vertexCount 'empty'             ==  0
-- vertexCount ('vertex' x)        ==  1
-- vertexCount                   ==  'length' . 'vertexList'
-- vertexCount x \< vertexCount y ==> x \< y
-- @
vertexCount :: AdjacencyIntMap e -> Int
vertexCount = M.size . adjacencyIntMap

-- | The number of (non-'zero') edges in a graph.
-- Complexity: /O(n)/ time.
--
-- @
-- edgeCount 'empty'        == 0
-- edgeCount ('vertex' x)   == 0
-- edgeCount ('edge' e x y) == if e == 'zero' then 0 else 1
-- edgeCount              == 'length' . 'edgeList'
-- @
edgeCount :: AdjacencyIntMap e -> Int
edgeCount = getSum . foldMap (Sum . M.size) . adjacencyIntMap

-- | The sorted list of vertices of a given graph.
-- Complexity: /O(n)/ time and memory.
--
-- @
-- vertexList 'empty'      == []
-- vertexList ('vertex' x) == [x]
-- vertexList . 'vertices' == 'Data.List.nub' . 'Data.List.sort'
-- @
vertexList :: AdjacencyIntMap e -> [Int]
vertexList = M.keys . adjacencyIntMap

-- | The list of edges of a graph, sorted lexicographically with respect to
-- pairs of connected vertices (i.e. edge-labels are ignored when sorting).
-- Complexity: /O(n + m)/ time and /O(m)/ memory.
--
-- @
-- edgeList 'empty'        == []
-- edgeList ('vertex' x)   == []
-- edgeList ('edge' e x y) == if e == 'zero' then [] else [(e,x,y)]
-- @
edgeList :: AdjacencyIntMap e -> [(e, Int, Int)]
edgeList (AM m) =
  [(e, x, y) | (x, ys) <- M.toAscList m, (y, e) <- M.toAscList ys]

-- | The set of vertices of a given graph.
-- Complexity: /O(n)/ time and memory.
--
-- @
-- vertexSet 'empty'      == Set.'Set.empty'
-- vertexSet . 'vertex'   == Set.'Set.singleton'
-- vertexSet . 'vertices' == Set.'Set.fromList'
-- @
vertexSet :: AdjacencyIntMap e -> IntSet
vertexSet = M.keysSet . adjacencyIntMap

-- | The set of edges of a given graph.
-- Complexity: /O(n + m)/ time and /O(m)/ memory.
--
-- @
-- edgeSet 'empty'        == Set.'Set.empty'
-- edgeSet ('vertex' x)   == Set.'Set.empty'
-- edgeSet ('edge' e x y) == if e == 'zero' then Set.'Set.empty' else Set.'Set.singleton' (e,x,y)
-- @
edgeSet :: (Eq e) => AdjacencyIntMap e -> Set (e, Int, Int)
edgeSet = S.fromAscList . edgeList

-- | The /preset/ of an element @x@ is the set of its /direct predecessors/.
-- Complexity: /O(n * log(n))/ time and /O(n)/ memory.
--
-- @
-- preSet x 'empty'        == Set.'Set.empty'
-- preSet x ('vertex' x)   == Set.'Set.empty'
-- preSet 1 ('edge' e 1 2) == Set.'Set.empty'
-- preSet y ('edge' e x y) == if e == 'zero' then Set.'Set.empty' else Set.'Set.fromList' [x]
-- @
preSet :: Int -> AdjacencyIntMap e -> IntSet
preSet x (AM m) =
  IS.fromAscList
    [a | (a, es) <- M.toAscList m, M.member x es]

-- | The /postset/ of a vertex is the set of its /direct successors/.
-- Complexity: /O(log(n))/ time and /O(1)/ memory.
--
-- @
-- postSet x 'empty'        == Set.'Set.empty'
-- postSet x ('vertex' x)   == Set.'Set.empty'
-- postSet x ('edge' e x y) == if e == 'zero' then Set.'Set.empty' else Set.'Set.fromList' [y]
-- postSet 2 ('edge' e 1 2) == Set.'Set.empty'
-- @
postSet :: Int -> AdjacencyIntMap e -> IntSet
postSet x = M.keysSet . M.findWithDefault M.empty x . adjacencyIntMap

-- TODO: Optimise.

-- | Convert a graph to the corresponding unlabelled 'AM.AdjacencyIntMap' by
-- forgetting labels on all non-'zero' edges.
-- Complexity: /O((n + m) * log(n))/ time and memory.
--
-- @
-- 'hasEdge' x y == 'AM.hasEdge' x y . skeleton
-- @
skeleton :: AdjacencyIntMap e -> AM.AdjacencyIntMap
skeleton (AM m) = AM.fromAdjacencyIntSets $ M.toAscList $ M.map M.keysSet m

-- | Remove a vertex from a given graph.
-- Complexity: /O(n*log(n))/ time.
--
-- @
-- removeVertex x ('vertex' x)       == 'empty'
-- removeVertex 1 ('vertex' 2)       == 'vertex' 2
-- removeVertex x ('edge' e x x)     == 'empty'
-- removeVertex 1 ('edge' e 1 2)     == 'vertex' 2
-- removeVertex x . removeVertex x == removeVertex x
-- @
removeVertex :: Int -> AdjacencyIntMap e -> AdjacencyIntMap e
removeVertex x = AM . M.map (M.delete x) . M.delete x . adjacencyIntMap

-- | Remove an edge from a given graph.
-- Complexity: /O(log(n))/ time.
--
-- @
-- removeEdge x y ('edge' e x y)     == 'vertices' [x,y]
-- removeEdge x y . removeEdge x y == removeEdge x y
-- removeEdge x y . 'removeVertex' x == 'removeVertex' x
-- removeEdge 1 1 (1 * 1 * 2 * 2)  == 1 * 2 * 2
-- removeEdge 1 2 (1 * 1 * 2 * 2)  == 1 * 1 + 2 * 2
-- @
removeEdge :: Int -> Int -> AdjacencyIntMap e -> AdjacencyIntMap e
removeEdge x y = AM . M.adjust (M.delete y) x . adjacencyIntMap

-- | The function @'replaceVertex' x y@ replaces vertex @x@ with vertex @y@ in a
-- given 'AdjacencyIntMap'. If @y@ already exists, @x@ and @y@ will be merged.
-- Complexity: /O((n + m) * log(n))/ time.
--
-- @
-- replaceVertex x x            == id
-- replaceVertex x y ('vertex' x) == 'vertex' y
-- replaceVertex x y            == 'gmap' (\\v -> if v == x then y else v)
-- @
replaceVertex :: (Eq e, Monoid e) => Int -> Int -> AdjacencyIntMap e -> AdjacencyIntMap e
replaceVertex u v = gmap $ \w -> if w == u then v else w

-- | Replace an edge from a given graph. If it doesn't exist, it will be created.
-- Complexity: /O(log(n))/ time.
--
-- @
-- replaceEdge e x y z                 == 'overlay' (removeEdge x y z) ('edge' e x y)
-- replaceEdge e x y ('edge' f x y)      == 'edge' e x y
-- 'edgeLabel' x y (replaceEdge e x y z) == e
-- @
replaceEdge :: (Eq e, Monoid e) => e -> Int -> Int -> AdjacencyIntMap e -> AdjacencyIntMap e
replaceEdge e x y
  | e == zero = AM . addY . M.alter (Just . maybe M.empty (M.delete y)) x . adjacencyIntMap
  | otherwise = AM . addY . M.alter replace x . adjacencyIntMap
  where
    addY = M.alter (Just . fromMaybe M.empty) y
    replace (Just m) = Just $ M.insert y e m
    replace Nothing = Just $ M.singleton y e

-- | Transpose a given graph.
-- Complexity: /O(m * log(n))/ time, /O(n + m)/ memory.
--
-- @
-- transpose 'empty'        == 'empty'
-- transpose ('vertex' x)   == 'vertex' x
-- transpose ('edge' e x y) == 'edge' e y x
-- transpose . transpose  == id
-- @
transpose :: Monoid e => AdjacencyIntMap e -> AdjacencyIntMap e
transpose (AM m) = AM $ M.foldrWithKey combine vs m
  where
    -- No need to use @nonZeroUnion@ here, since we do not add any new edges
    combine v es =
      M.unionWith (M.unionWith mappend) $
        M.fromAscList [(u, M.singleton v e) | (u, e) <- M.toAscList es]
    vs = M.fromSet (const M.empty) (M.keysSet m)

-- | Transform a graph by applying a function to each of its vertices. This is
-- similar to @Functor@'s 'fmap' but can be used with non-fully-parametric
-- 'AdjacencyIntMap'.
-- Complexity: /O((n + m) * log(n))/ time.
--
-- @
-- gmap f 'empty'        == 'empty'
-- gmap f ('vertex' x)   == 'vertex' (f x)
-- gmap f ('edge' e x y) == 'edge' e (f x) (f y)
-- gmap 'id'             == 'id'
-- gmap f . gmap g     == gmap (f . g)
-- @
gmap :: (Eq e, Monoid e) => (Int -> Int) -> AdjacencyIntMap e -> AdjacencyIntMap e
gmap f =
  AM . trimZeroes . M.map (M.mapKeysWith mappend f)
    . M.mapKeysWith (M.unionWith mappend) f
    . adjacencyIntMap

-- | Transform a graph by applying a function @h@ to each of its edge labels.
-- Complexity: /O((n + m) * log(n))/ time.
--
-- The function @h@ is required to be a /homomorphism/ on the underlying type of
-- labels @e@. At the very least it must preserve 'zero' and '<+>':
--
-- @
-- h 'zero'      == 'zero'
-- h x '<+>' h y == h (x '<+>' y)
-- @
--
-- If @e@ is also a semiring, then @h@ must also preserve the multiplicative
-- structure:
--
-- @
-- h 'one'       == 'one'
-- h x '<.>' h y == h (x '<.>' y)
-- @
--
-- If the above requirements hold, then the implementation provides the
-- following guarantees.
--
-- @
-- emap h 'empty'           == 'empty'
-- emap h ('vertex' x)      == 'vertex' x
-- emap h ('edge' e x y)    == 'edge' (h e) x y
-- emap h ('overlay' x y)   == 'overlay' (emap h x) (emap h y)
-- emap h ('connect' e x y) == 'connect' (h e) (emap h x) (emap h y)
-- emap 'id'                == 'id'
-- emap g . emap h        == emap (g . h)
-- @
emap :: (Eq f, Monoid f) => (e -> f) -> AdjacencyIntMap e -> AdjacencyIntMap f
emap h = AM . trimZeroes . M.map (M.map h) . adjacencyIntMap

-- | Construct the /induced subgraph/ of a given graph by removing the
-- vertices that do not satisfy a given predicate.
-- Complexity: /O(n + m)/ time, assuming that the predicate takes /O(1)/ to
-- be evaluated.
--
-- @
-- induce ('const' True ) x      == x
-- induce ('const' False) x      == 'empty'
-- induce (/= x)               == 'removeVertex' x
-- induce p . induce q         == induce (\\x -> p x && q x)
-- 'isSubgraphOf' (induce p x) x == True
-- @
induce :: (Int -> Bool) -> AdjacencyIntMap e -> AdjacencyIntMap e
induce p =
  AM . M.map (M.filterWithKey (\k _ -> p k))
    . M.filterWithKey (\k _ -> p k)
    . adjacencyIntMap

-- -- | Construct the /induced subgraph/ of a given graph by removing the vertices
-- -- that are 'Nothing'.
-- -- Complexity: /O(n + m)/ time.
-- --
-- -- @
-- -- induceJust ('vertex' 'Nothing')                               == 'empty'
-- -- induceJust ('edge' ('Just' x) 'Nothing')                        == 'vertex' x
-- -- induceJust . 'gmap' 'Just'                                    == 'id'
-- -- induceJust . 'gmap' (\\x -> if p x then 'Just' x else 'Nothing') == 'induce' p
-- -- @
-- induceJust :: Ord a => AdjacencyIntMap e (Maybe a) -> AdjacencyIntMap e
-- induceJust = AM . M.map catMaybesMap . catMaybesMap . adjacencyIntMap
--   where
--     catMaybesMap = M.mapKeysMonotonic fromJust . M.delete Nothing

-- | Compute the /reflexive and transitive closure/ of a graph over the
-- underlying star semiring using the Warshall-Floyd-Kleene algorithm.
--
-- @
-- closure 'empty'         == 'empty'
-- closure ('vertex' x)    == 'edge' 'one' x x
-- closure ('edge' e x x)  == 'edge' 'one' x x
-- closure ('edge' e x y)  == 'edges' [('one',x,x), (e,x,y), ('one',y,y)]
-- closure               == 'reflexiveClosure' . 'transitiveClosure'
-- closure               == 'transitiveClosure' . 'reflexiveClosure'
-- closure . closure     == closure
-- 'postSet' x (closure y) == Set.'Set.fromList' ('Algebra.Graph.ToGraph.reachable' x y)
-- @
closure :: (Eq e, StarSemiring e) => AdjacencyIntMap e -> AdjacencyIntMap e
closure = goWarshallFloydKleene . reflexiveClosure

-- | Compute the /reflexive closure/ of a graph over the underlying semiring by
-- adding a self-loop of weight 'one' to every vertex.
-- Complexity: /O(n * log(n))/ time.
--
-- @
-- reflexiveClosure 'empty'              == 'empty'
-- reflexiveClosure ('vertex' x)         == 'edge' 'one' x x
-- reflexiveClosure ('edge' e x x)       == 'edge' 'one' x x
-- reflexiveClosure ('edge' e x y)       == 'edges' [('one',x,x), (e,x,y), ('one',y,y)]
-- reflexiveClosure . reflexiveClosure == reflexiveClosure
-- @
reflexiveClosure :: (Semiring e) => AdjacencyIntMap e -> AdjacencyIntMap e
reflexiveClosure (AM m) = AM $ M.mapWithKey (\k -> M.insertWith (<+>) k one) m

-- | Compute the /symmetric closure/ of a graph by overlaying it with its own
-- transpose.
-- Complexity: /O((n + m) * log(n))/ time.
--
-- @
-- symmetricClosure 'empty'              == 'empty'
-- symmetricClosure ('vertex' x)         == 'vertex' x
-- symmetricClosure ('edge' e x y)       == 'edges' [(e,x,y), (e,y,x)]
-- symmetricClosure x                  == 'overlay' x ('transpose' x)
-- symmetricClosure . symmetricClosure == symmetricClosure
-- @
symmetricClosure :: (Eq e, Monoid e) => AdjacencyIntMap e -> AdjacencyIntMap e
symmetricClosure m = overlay m (transpose m)

-- | Compute the /transitive closure/ of a graph over the underlying star
-- semiring using a modified version of the Warshall-Floyd-Kleene algorithm,
-- which omits the reflexivity step.
--
-- @
-- transitiveClosure 'empty'               == 'empty'
-- transitiveClosure ('vertex' x)          == 'vertex' x
-- transitiveClosure ('edge' e x y)        == 'edge' e x y
-- transitiveClosure . transitiveClosure == transitiveClosure
-- @
transitiveClosure :: (Eq e, StarSemiring e) => AdjacencyIntMap e -> AdjacencyIntMap e
transitiveClosure = goWarshallFloydKleene

-- The iterative part of the Warshall-Floyd-Kleene algorithm
goWarshallFloydKleene :: (Eq e, StarSemiring e) => AdjacencyIntMap e -> AdjacencyIntMap e
goWarshallFloydKleene (AM m) = AM $ foldr update m vs
  where
    vs = IS.toAscList (M.keysSet m)
    update k cur = M.fromAscList [(i, go i (get i k <.> starkk)) | i <- vs]
      where
        get i j = edgeLabel i j (AM cur)
        starkk = star (get k k)
        go i ik =
          M.fromAscList
            [(j, e) | j <- vs, let e = get i j <+> ik <.> get k j, e /= zero]

-- | Check that the internal graph representation is consistent, i.e. that all
-- edges refer to existing vertices, and there are no 'zero'-labelled edges. It
-- should be impossible to create an inconsistent adjacency map, and we use this
-- function in testing.
consistent :: (Eq e, Monoid e) => AdjacencyIntMap e -> Bool
consistent (AM m) =
  referredToVertexSet m `IS.isSubsetOf` M.keysSet m
    && and [e /= zero | (_, es) <- M.toAscList m, (_, e) <- M.toAscList es]

-- The set of vertices that are referred to by the edges in an adjacency map
referredToVertexSet :: IntMap (IntMap e) -> IntSet
referredToVertexSet m =
  IS.fromList $
    concat
      [[x, y] | (x, ys) <- M.toAscList m, (y, _) <- M.toAscList ys]
