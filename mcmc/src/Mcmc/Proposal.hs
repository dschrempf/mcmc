{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes #-}

-- TODO: Proposals on simplices: SimplexElementScale (?).

-- TODO: Proposals on trees:
-- - Slide a node on the tree.
-- - Scale a tree.

-- TODO: Proposals on tree topologies.
-- - NNI
-- - Narrow (what is this, see RevBayes)
-- - FNPR (dito)

-- |
-- Module      :  Mcmc.Proposal
-- Description :  Proposals and cycles
-- Copyright   :  (c) Dominik Schrempf 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Wed May 20 13:42:53 2020.
module Mcmc.Proposal
  ( -- * Proposal
    Proposal (..),
    (>>>),
    ProposalSimple (..),
    Tuner (tParam, tFunc),
    tuner,
    tune,

    -- * Cycle
    Order (..),
    Cycle (ccProposals),
    fromList,
    setOrder,
    getNCycles,
    tuneCycle,
    autotuneCycle,
    summarizeCycle,

    -- * Acceptance
    Acceptance (fromAcceptance),
    emptyA,
    pushA,
    resetA,
    transformKeysA,
    acceptanceRatios,
  )
where

import Data.Aeson
import Data.Default
import Data.Function
import Data.List
import qualified Data.Map.Strict as M
import Data.Map.Strict (Map)
import Data.Maybe
import qualified Data.Text.Lazy as T
import Data.Text.Lazy (Text)
import qualified Data.Text.Lazy.Builder as B
import qualified Data.Text.Lazy.Builder.Int as B
import qualified Data.Text.Lazy.Builder.RealFloat as B
import Lens.Micro
import Mcmc.Tools.Shuffle
import Numeric.Log hiding (sum)
import System.Random.MWC

-- | A 'Proposal' is an instruction about how the Markov chain will traverse the
-- state space @a@. Essentially, it is a probability mass or probability density
-- conditioned on the current state (i.e., a kernel).
--
-- A 'Proposal' may be tuneable in that it contains information about how to enlarge
-- or shrink the step size to tune the acceptance ratio.
data Proposal a
  = Proposal
      { -- | Name (no proposals with the same name are allowed in a 'Cycle').
        pName :: String,
        -- | The weight determines how often a 'Proposal' is executed per iteration of
        -- the Markov chain.
        pWeight :: Int,
        -- | Simple proposal without tuning information.
        pSimple :: ProposalSimple a,
        -- | Tuning is disabled if set to 'Nothing'.
        pTuner :: Maybe (Tuner a)
      }

instance Show (Proposal a) where
  show m = show $ pName m

instance Eq (Proposal a) where
  m == n = pName m == pName n

instance Ord (Proposal a) where
  compare = compare `on` pName

convertP :: Lens' b a -> Proposal a -> Proposal b
convertP l (Proposal n w s t) = Proposal n w (convertS l s) (convertT l <$> t)

-- | Convert a proposal from one data type to another using a lens.
--
-- For example:
--
-- @
-- scaleFirstEntryOfTuple = scale >>> _1
-- @
(>>>) :: Lens' b a -> Proposal a -> Proposal b
(>>>) = convertP

-- One could also use a different type for 'pSample', so that 'pKernel' can
-- be avoided. In detail,
--
-- @
--   pSample :: a -> GenIO -> IO (a, Log Double, Log, Double)
-- @
--
-- where the kernels describe the probability of going there and back. However,
-- we may need more information about the proposal for other MCMC samplers
-- different from Metropolis-Hastings.

-- | Simple proposal without tuning information.
--
-- In order to calculate the Metropolis-Hastings ratio, we need to know the
-- kernel (i.e., the probability mass or probability density) of jumping
-- forwards and backwards.
data ProposalSimple a
  = ProposalSimple
      { -- | Instruction about randomly moving from the current state to a new
        -- state, given some source of randomness.
        pSample :: a -> GenIO -> IO a,
        -- | The kernel of going from one state to another. Set to 'Nothing' for
        -- symmetric proposals.
        pKernel :: Maybe (a -> a -> Log Double)
      }

convertS :: Lens' b a -> ProposalSimple a -> ProposalSimple b
convertS l (ProposalSimple s mk) = ProposalSimple s' mk'
  where s' v g = do x' <- s (v ^. l) g
                    return $ set l x' v
        mk' = case mk of
          Nothing -> Nothing
          Just k -> Just $ \x y -> k (x ^. l) (y ^. l)

-- | Tune the acceptance ratio of a 'Proposal'; see 'tune', or 'autotuneCycle'.
data Tuner a
  = Tuner
      { tParam :: Double,
        tFunc :: Double -> ProposalSimple a
      }

convertT :: Lens' b a -> Tuner a -> Tuner b
convertT l (Tuner p f) = Tuner p f'
  where f' x = convertS l $ f x

-- | Create a 'Tuner'. The tuning function accepts a tuning parameter, and
-- returns a corresponding 'ProposalSimple'. The larger the tuning parameter, the
-- larger the 'Proposal', and vice versa.
tuner :: (Double -> ProposalSimple a) -> Tuner a
tuner = Tuner 1.0

-- Minimal tuning parameter; subject to change.
tuningParamMin :: Double
tuningParamMin = 1e-12

-- | Tune a 'Proposal'. Return 'Nothing' if 'Proposal' is not tuneable. If the parameter
--   @dt@ is larger than 1.0, the 'Proposal' is enlarged, if @0<dt<1.0@, it is
--   shrunk. Negative tuning parameters are not allowed.
tune :: Double -> Proposal a -> Maybe (Proposal a)
tune dt m
  | dt <= 0 = error $ "tune: Tuning parameter not positive: " <> show dt <> "."
  | otherwise = do
    (Tuner t f) <- pTuner m
    -- Ensure that the tuning parameter is not too small.
    let t' = max tuningParamMin (t * dt)
    return $ m {pSimple = f t', pTuner = Just $ Tuner t' f}

-- XXX: The desired acceptance ratio 0.44 is optimal for one-dimensional
-- 'Proposal's; one could also store the affected number of dimensions with the
-- 'Proposal' and tune towards an acceptance ratio accounting for the number of
-- dimensions.
ratioOpt :: Double
ratioOpt = 0.44

-- | Define the order in which 'Proposal's are executed in a 'Cycle'. The total
-- number of 'Proposal's per 'Cycle' may differ between 'Order's (e.g., compare
-- 'RandomO' and 'RandomReversibleO').
data Order
  = -- | Shuffle the 'Proposal's in the 'Cycle'. The 'Proposal's are replicated
    -- according to their weights and executed in random order. If a 'Proposal' has
    -- weight @w@, it is executed exactly @w@ times per iteration.
    RandomO
  | -- | The 'Proposal's are executed sequentially, in the order they appear in the
    -- 'Cycle'. 'Proposal's with weight @w>1@ are repeated immediately @w@ times
    -- (and not appended to the end of the list).
    SequentialO
  | -- | Similar to 'RandomO'. However, a reversed copy of the list of
    --  shuffled 'Proposal's is appended such that the resulting Markov chain is
    --  reversible.
    --  Note: the total number of 'Proposal's executed per cycle is twice the number
    --  of 'RandomO'.
    RandomReversibleO
  | -- | Similar to 'SequentialO'. However, a reversed copy of the list of
    -- sequentially ordered 'Proposal's is appended such that the resulting Markov
    -- chain is reversible.
    SequentialReversibleO
  deriving (Eq, Show)

instance Default Order where def = RandomO

-- | In brief, a 'Cycle' is a list of proposals. The state of the Markov chain will
-- be logged only after all 'Proposal's in the 'Cycle' have been completed, and the
-- iteration counter will be increased by one. The order in which the 'Proposal's
-- are executed is specified by 'Order'. The default is 'RandomO'.
--
-- __Proposals must have unique names__, so that they can be identified.
data Cycle a
  = Cycle
      { ccProposals :: [Proposal a],
        ccOrder :: Order
      }

-- | Create a 'Cycle' from a list of 'Proposal's.
fromList :: [Proposal a] -> Cycle a
fromList [] =
  error "fromList: Received an empty list but cannot create an empty Cycle."
fromList xs =
  if length (nub nms) == length nms
    then Cycle xs def
    else error "fromList: Proposals don't have unique names."
  where
    nms = map pName xs

-- | Set the order of 'Proposal's in a 'Cycle'.
setOrder :: Order -> Cycle a -> Cycle a
setOrder o c = c {ccOrder = o}

-- | Replicate 'Proposal's according to their weights and possibly shuffle them.
getNCycles :: Cycle a -> Int -> GenIO -> IO [[Proposal a]]
getNCycles (Cycle xs o) n g = case o of
  RandomO -> shuffleN ps n g
  SequentialO -> return $ replicate n ps
  RandomReversibleO -> do
    psRs <- shuffleN ps n g
    return [psR ++ reverse psR | psR <- psRs]
  SequentialReversibleO -> return $ replicate n $ ps ++ reverse ps
  where
    !ps = concat [replicate (pWeight m) m | m <- xs]

-- | Tune 'Proposal's in the 'Cycle'. See 'tune'.
tuneCycle :: Map (Proposal a) Double -> Cycle a -> Cycle a
tuneCycle m c =
  if sort (M.keys m) == sort ps
    then c {ccProposals = map tuneF ps}
    else error "tuneCycle: Map contains proposals that are not in the cycle."
  where
    ps = ccProposals c
    tuneF p = case m M.!? p of
      Nothing -> p
      Just x -> fromMaybe p (tune x p)

-- | Calculate acceptance ratios and auto tune the 'Proposal's in the 'Cycle'. For
-- now, a 'Proposal' is enlarged when the acceptance ratio is above 0.44, and
-- shrunk otherwise. Do not change 'Proposal's that are not tuneable.
autotuneCycle :: Acceptance (Proposal a) -> Cycle a -> Cycle a
autotuneCycle a = tuneCycle (M.map (\x -> exp $ x - ratioOpt) $ acceptanceRatios a)

renderRow :: Text -> Text -> Text -> Text -> Text -> Text -> Text
renderRow name weight nAccept nReject acceptRatio tuneParam = "   " <> nm <> wt <> na <> nr <> ra <> tp
  where
    nm = T.justifyLeft 30 ' ' name
    wt = T.justifyRight 8 ' ' weight
    na = T.justifyRight 15 ' ' nAccept
    nr = T.justifyRight 15 ' ' nReject
    ra = T.justifyRight 15 ' ' acceptRatio
    tp = T.justifyRight 20 ' ' tuneParam

proposalHeader :: Text
proposalHeader =
  renderRow "Proposal" "Weight" "Accepted" "Rejected" "Ratio" "Tuning parameter"

summarizeProposal :: Proposal a -> Maybe (Int, Int, Double) -> Text
summarizeProposal m r = renderRow (T.pack name) weight nAccept nReject acceptRatio tuneParamStr
  where
    name = pName m
    weight = B.toLazyText $ B.decimal $ pWeight m
    nAccept = B.toLazyText $ maybe "" (B.decimal . (^. _1)) r
    nReject = B.toLazyText $ maybe "" (B.decimal . (^. _2)) r
    acceptRatio = B.toLazyText $ maybe "" (B.formatRealFloat B.Fixed (Just 3) . (^. _3)) r
    tuneParamStr = B.toLazyText $ maybe "" (B.formatRealFloat B.Fixed (Just 3)) (tParam <$> pTuner m)

-- | Summarize the 'Proposal's in the 'Cycle'. Also report acceptance ratios.
summarizeCycle :: Acceptance (Proposal a) -> Cycle a -> Text
summarizeCycle a c =
  T.intercalate "\n" $
    [ "Summary of proposal(s) in cycle. " <> mpi <> " proposal(s) per iteration.",
      proposalHeader,
      "   " <> T.replicate (T.length proposalHeader - 3) "─"
    ]
      ++ [summarizeProposal m (ar m) | m <- ps]
      ++ ["   " <> T.replicate (T.length proposalHeader - 3) "─"]
  where
    ps = ccProposals c
    mpi = B.toLazyText $ B.decimal $ sum $ map pWeight ps
    ar m = acceptanceRatio m a

-- | For each key @k@, store the number of accepted and rejected proposals.
newtype Acceptance k = Acceptance {fromAcceptance :: Map k (Int, Int)}

instance ToJSONKey k => ToJSON (Acceptance k) where
  toJSON (Acceptance m) = toJSON m
  toEncoding (Acceptance m) = toEncoding m

instance (Ord k, FromJSONKey k) => FromJSON (Acceptance k) where
  parseJSON v = Acceptance <$> parseJSON v

-- | In the beginning there was the Word.
--
-- Initialize an empty storage of accepted/rejected values.
emptyA :: Ord k => [k] -> Acceptance k
emptyA ks = Acceptance $ M.fromList [(k, (0, 0)) | k <- ks]

-- | For key @k@, prepend an accepted (True) or rejected (False) proposal.
pushA :: (Ord k, Show k) => k -> Bool -> Acceptance k -> Acceptance k
pushA k True = Acceptance . M.adjust (\(a, r) -> (succ a, r)) k . fromAcceptance
pushA k False = Acceptance . M.adjust (\(a, r) -> (a, succ r)) k . fromAcceptance
{-# INLINEABLE pushA #-}

-- | Reset acceptance storage.
resetA :: Ord k => Acceptance k -> Acceptance k
resetA = emptyA . M.keys . fromAcceptance

transformKeys :: (Ord k1, Ord k2) => [k1] -> [k2] -> Map k1 v -> Map k2 v
transformKeys ks1 ks2 m = foldl' insrt M.empty $ zip ks1 ks2
  where
    insrt m' (k1, k2) = M.insert k2 (m M.! k1) m'

-- | Transform keys using the given lists. Keys not provided will not be present
-- in the new 'Acceptance' variable.
transformKeysA :: (Ord k1, Ord k2) => [k1] -> [k2] -> Acceptance k1 -> Acceptance k2
transformKeysA ks1 ks2 = Acceptance . transformKeys ks1 ks2 . fromAcceptance

-- | Acceptance counts and ratio for a specific proposal.
acceptanceRatio :: (Show k, Ord k) => k -> Acceptance k -> Maybe (Int, Int, Double)
acceptanceRatio k a = case fromAcceptance a M.!? k of
  Just (0, 0) -> Nothing
  Just (as, rs) -> Just (as, rs, fromIntegral as / fromIntegral (as + rs))
  Nothing -> error $ "acceptanceRatio: Key not found in map: " ++ show k ++ "."

-- | Acceptance ratios for all proposals.
acceptanceRatios :: Acceptance k -> Map k Double
acceptanceRatios = M.map (\(as, rs) -> fromIntegral as / fromIntegral (as + rs)) . fromAcceptance
