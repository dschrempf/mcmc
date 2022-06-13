-- |
-- Module      :  Mcmc.Proposal.Hamiltonian.Internal
-- Description :  Internal definitions related to Hamiltonian dynamics
-- Copyright   :  2022 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Thu Jun  9 15:12:39 2022.
--
-- See "Mcmc.Proposal.Hamiltonian.Hamiltonian".
--
-- References:
--
-- - [1] Chapter 5 of Handbook of Monte Carlo: Neal, R. M., MCMC Using
--   Hamiltonian Dynamics, In S. Brooks, A. Gelman, G. Jones, & X. Meng (Eds.),
--   Handbook of Markov Chain Monte Carlo (2011), CRC press.
--
-- - [2] Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B., Bayesian data
--   analysis (2014), CRC Press.
--
-- - [3] Review by Betancourt and notes: Betancourt, M., A conceptual
--   introduction to Hamiltonian Monte Carlo, arXiv, 1701–02434 (2017).
--
-- - [4] Matthew D. Hoffman, Andrew Gelman (2014) The No-U-Turn Sampler:
--   Adaptively Setting Path Lengths in Hamiltonian Monte Carlo, Journal of
--   Machine Learning Research.
module Mcmc.Proposal.Hamiltonian.Internal
  ( -- * Parameters
    Mu,
    MassesInv,
    HData (..),
    getHData,
    HParamsI (..),
    hParamsIWith,

    -- * Tuning
    Dimension,
    massesToVector,
    vectorToMasses,
    toAuxiliaryTuningParameters,
    fromAuxiliaryTuningParameters,
    findReasonableEpsilon,
    hTuningFunctionWith,

    -- * Structure of state
    checkHStructureWith,

    -- * Hamiltonian dynamics
    generateMomenta,
    exponentialKineticEnergy,

    -- * Leapfrog integrator
    Target,
    leapfrog,
  )
where

import Control.Monad
import Control.Monad.ST
import Data.Foldable
import qualified Data.Vector as VB
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as VU
import Mcmc.Proposal
import Mcmc.Proposal.Hamiltonian.Common
import qualified Numeric.LinearAlgebra as L
import Numeric.Log
import qualified Statistics.Covariance as S
import qualified Statistics.Function as S
import qualified Statistics.Sample as S
import System.Random.MWC
import System.Random.Stateful

-- Variable tuning parameters.
--
-- See Algorithm 5 or 6 in [4].
data TParamsVar = TParamsVar
  { -- \bar{eps} of Algorithm 5 or 6.
    tpvLeapfrogScalingFactorMean :: LeapfrogScalingFactor,
    -- H_i of Algorithm 5 or 6.
    tpvHStatistics :: Double,
    -- m of Algorithm 5 or 6.
    tpvCurrentTuningStep :: Double
  }
  deriving (Show)

tParamsVar :: TParamsVar
tParamsVar = TParamsVar 1.0 0.0 1.0

-- Fixed tuning parameters.
--
-- See Algorithm 5 and 6 in [4].
data TParamsFixed = TParamsFixed
  { tpfEps0 :: Double,
    tpfMu :: Double,
    tpfGa :: Double,
    tpfT0 :: Double,
    tpfKa :: Double
  }
  deriving (Show)

-- The default tuning parameters in [4] are:
--
--   mu = log $ 10 * eps
--   ga = 0.05
--   t0 = 10
--   ka = 0.75
--
-- However, these default tuning parameters will be off, because the authors
-- suggesting these values tune the proposal after every single iteration.
--
-- The following values are tweaked for our case, where tuning does not happen
-- after each iteration. Of course, we could tune the leapfrog parameters after
-- each generation. Even the mass parameters could be tuned each iteration when
-- the masses are estimated from more past iterations spanning many tuning
-- intervals.
--
-- NOTE: In theory, these we could expose these internal tuning parameters to
-- the user.
tParamsFixedWith :: LeapfrogScalingFactor -> TParamsFixed
tParamsFixedWith eps = TParamsFixed eps mu ga t0 ka
  where
    mu = log $ 10 * eps
    ga = 0.1
    t0 = 3
    ka = 0.5

-- Mean vector containing zeroes. We save this vector because it is required
-- when sampling from the multivariate normal distribution.
type Mu = L.Vector Double

-- Symmetric, inverted mass matrix.
type MassesInv = L.Herm Double

-- Data type containing memoized values.
data HData = HData
  { hMu :: Mu,
    hMassesInv :: MassesInv
  }
  deriving (Show)

-- Compute inverted mass matrix.
--
-- Call 'error' if the determinant of the covariance matrix is negative.
getHData :: Masses -> HData
getHData ms =
  -- The multivariate normal distribution requires a positive definite matrix
  -- with positive determinant.
  if sign == 1.0
    then HData mu massesInvH
    else
      let msg =
            "getHData: Determinant of covariance matrix is negative."
       in error msg
  where
    nrows = L.rows $ L.unSym ms
    mu = L.fromList $ replicate nrows 0.0
    (massesInv, (_, sign)) = L.invlndet $ L.unSym ms
    -- NOTE: In theory we can trust that the matrix is symmetric here, because
    -- the inverse of a symmetric matrix is symmetric. However, one may want to
    -- implement a check anyways.
    massesInvH = L.trustSym massesInv

-- All internal parameters.
data HParamsI = HParamsI
  { hpsLeapfrogScalingFactor :: LeapfrogScalingFactor,
    hpsLeapfrogSimulationLength :: LeapfrogSimulationLength,
    hpsMasses :: Masses,
    hpsTParamsVar :: TParamsVar,
    hpsTParamsFixed :: TParamsFixed,
    hpsData :: HData
  }
  deriving (Show)

-- NOTE: If changed, amend help text of 'defaultHParams', and 'defaultNParams'.
defaultLeapfrogScalingFactor :: LeapfrogScalingFactor
defaultLeapfrogScalingFactor = 0.1

-- NOTE: If changed, amend help text of 'defaultHParams'.
defaultLeapfrogSimulationLength :: LeapfrogSimulationLength
defaultLeapfrogSimulationLength = 0.5

-- Instantiate all internal parameters.
hParamsIWith ::
  Target ->
  Positions ->
  Maybe LeapfrogScalingFactor ->
  Maybe LeapfrogSimulationLength ->
  Maybe Masses ->
  Either String HParamsI
hParamsIWith htarget p mEps mLa mMs = do
  d <- case VS.length p of
    0 -> eWith "Empty position vector."
    d -> Right d
  let defMs = L.trustSym $ L.ident d
  ms <- case mMs of
    Nothing -> Right defMs
    Just ms -> do
      let ms' = L.unSym ms
          diagonalMs = L.toList $ L.takeDiag ms'
      when (any (<= 0) diagonalMs) $ eWith "Some diagonal masses are zero or negative."
      let nrows = L.rows ms'
          ncols = L.cols ms'
      when (nrows /= ncols) $ eWith "Mass matrix is not square."
      Right ms
  let hdata = getHData ms
  la <- case mLa of
    Nothing -> Right defaultLeapfrogSimulationLength
    Just l
      | l <= 0 -> eWith "Leapfrog simulation length is zero or negative."
      | otherwise -> Right l
  eps <- case mEps of
    Nothing -> Right $ runST $ do
      g <- create
      findReasonableEpsilon htarget ms p g
    Just e
      | e <= 0 -> eWith "Leapfrog scaling factor is zero or negative."
      | otherwise -> Right e
  let tParamsFixed = tParamsFixedWith eps
  pure $ HParamsI eps la ms tParamsVar tParamsFixed hdata
  where
    eWith m = Left $ "hParamsIWith: " <> m

-- Dimension of the proposal.
type Dimension = Int

massesToVector :: Masses -> VU.Vector Double
massesToVector = VU.convert . L.flatten . L.unSym

vectorToMasses :: Dimension -> VU.Vector Double -> Masses
vectorToMasses d = L.trustSym . L.reshape d . VU.convert

-- Save internal parameters.
toAuxiliaryTuningParameters :: HParamsI -> AuxiliaryTuningParameters
toAuxiliaryTuningParameters (HParamsI eps la ms tpv tpf _) =
  -- Put masses to the end. Like so, conversion is easier.
  VU.fromList $ eps : la : epsMean : h : m : eps0 : mu : ga : t0 : ka : msL
  where
    (TParamsVar epsMean h m) = tpv
    (TParamsFixed eps0 mu ga t0 ka) = tpf
    msL = VU.toList $ massesToVector ms

-- Load internal parameters.
fromAuxiliaryTuningParameters :: Dimension -> AuxiliaryTuningParameters -> Either String HParamsI
fromAuxiliaryTuningParameters d xs
  | (d * d) + 10 /= len = Left "fromAuxiliaryTuningParameters: Dimension mismatch."
  | fromIntegral (d * d) /= lenMs = Left "fromAuxiliaryTuningParameters: Masses dimension mismatch."
  | otherwise = case VU.toList $ VU.take 10 xs of
      [eps, la, epsMean, h, m, eps0, mu, ga, t0, ka] ->
        let tpv = TParamsVar epsMean h m
            tpf = TParamsFixed eps0 mu ga t0 ka
         in Right $ HParamsI eps la ms tpv tpf hdata
      -- To please the exhaustive pattern match checker.
      _ -> Left "fromAuxiliaryTuningParameters: Impossible dimension mismatch."
  where
    len = VU.length xs
    msV = VU.drop 10 xs
    lenMs = VU.length msV
    ms = vectorToMasses d msV
    hdata = getHData ms

-- See Algorithm 4 in [4].
findReasonableEpsilon ::
  StatefulGen g m =>
  Target ->
  Masses ->
  Positions ->
  g ->
  m LeapfrogScalingFactor
findReasonableEpsilon t ms q g = do
  p <- generateMomenta mu ms g
  case leapfrog t msInv 1 eI q p of
    Nothing -> pure defaultLeapfrogScalingFactor
    Just (_, p', prQ, prQ') -> do
      let expEKin = exponentialKineticEnergy msInv p
          expEKin' = exponentialKineticEnergy msInv p'
          rI :: Double
          rI = exp $ ln $ prQ' * expEKin' / (prQ * expEKin)
          a :: Double
          a = if rI > 0.5 then 1 else (-1)
          go e r =
            if r ** a > 2 ** (negate a)
              then case leapfrog t msInv 1 e q p of
                Nothing -> e
                Just (_, p'', _, prQ'') ->
                  let expEKin'' = exponentialKineticEnergy msInv p''
                      r' :: Double
                      r' = exp $ ln $ prQ'' * expEKin'' / (prQ * expEKin)
                      e' = (2 ** a) * e
                   in go e' r'
              else e
      pure $ go eI rI
  where
    eI = 1.0
    (HData mu msInv) = getHData ms

-- If changed, also change help text of 'HTuneMasses'.
massMin :: Double
massMin = 1e-6

-- If changed, also change help text of 'HTuneMasses'.
massMax :: Double
massMax = 1e6

-- Minimal number of unique samples required for tuning the diagonal entries of
-- the mass matrix.
--
-- If changed, also change help text of 'HTuneMasses'.
samplesMinDiagonal :: Int
samplesMinDiagonal = 61

-- Minimal number of samples required for tuning all entries of the mass matrix.
--
-- If changed, also change help text of 'HTuneMasses'.
samplesMinAll :: Int
samplesMinAll = 201

getSampleSize :: VS.Vector Double -> Int
getSampleSize = VS.length . VS.uniq . S.gsort

-- Diagonal elements are variances which are strictly positive.
getNewMassDiagonalWithRescue :: Int -> Double -> Double -> Double
getNewMassDiagonalWithRescue sampleSize massOld massEstimate
  | sampleSize < samplesMinDiagonal = massOld
  -- NaN and negative masses could be errors.
  | isNaN massEstimate = massOld
  | massEstimate <= 0 = massOld
  | massMin > massNew = massMin
  | massNew > massMax = massMax
  | otherwise = massNew
  where
    massNewSqrt = recip 3 * (sqrt massOld + 2 * sqrt massEstimate)
    massNew = massNewSqrt ** 2

tuneDiagonalMassesOnly ::
  -- Conversion from value to vector.
  (a -> Positions) ->
  -- Value vector.
  VB.Vector a ->
  -- Old masses.
  Masses ->
  -- New masses.
  Masses
-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneDiagonalMassesOnly toVec xs ms
  -- If not enough data is available, do not tune.
  | VB.length xs < samplesMinDiagonal = ms
  | dimState /= dimMs = error "tuneDiagonalMassesOnly: Dimension mismatch."
  -- Replace the diagonal.
  | otherwise = L.trustSym $ msOld - L.diag msDiagonalOld + L.diag msDiagonalNew
  where
    -- xs: Each element contains all parameters of one iteration.
    -- xs': Each element is a vector containing all parameters changed by the
    -- proposal of one iteration.
    xs' = VB.map toVec xs
    -- xs'': Matrix with each row containing all parameter values changed by the
    -- proposal of one iteration.
    xs'' = L.fromRows $ VB.toList xs'
    -- We can safely use 'VB.head' here since the length of 'xs' must be larger
    -- than 'samplesMinDiagonal'.
    dimState = VS.length $ VB.head xs'
    sampleSizes = VS.fromList $ map getSampleSize $ L.toColumns xs''
    msOld = L.unSym ms
    dimMs = L.rows msOld
    msDiagonalOld = L.takeDiag msOld
    msDiagonalEstimate = VS.fromList $ map (recip . S.variance) $ L.toColumns xs''
    msDiagonalNew =
      VS.zipWith3
        getNewMassDiagonalWithRescue
        sampleSizes
        msDiagonalOld
        msDiagonalEstimate

tuneAllMasses ::
  -- Conversion from value to vector.
  (a -> Positions) ->
  -- Value vector.
  VB.Vector a ->
  -- Old masses.
  Masses ->
  -- New masses.
  Masses
-- NOTE: Here, we lose time because we convert the states to vectors again,
-- something that has already been done. But then, auto tuning is not a runtime
-- determining factor.
tuneAllMasses toVec xs ms
  -- If not enough data is available, do not tune.
  | VB.length xs < samplesMinDiagonal = ms
  -- If not enough data is available, only the diagonal masses are tuned.
  | VB.length xs < samplesMinAll = fallbackDiagonal
  | L.rank xs'' /= dimState = fallbackDiagonal
  | dimState /= dimMs = error "tuneAllMasses: Dimension mismatch."
  | otherwise = L.trustSym msNew
  where
    fallbackDiagonal = tuneDiagonalMassesOnly toVec xs ms
    -- xs: Each element contains all parameters of one iteration.
    -- xs': Each element is a vector containing all parameters changed by the
    -- proposal of one iteration.
    xs' = VB.map toVec xs
    -- xs'': Matrix with each row containing all parameter values changed by the
    -- proposal of one iteration.
    xs'' = L.fromRows $ VB.toList xs'
    -- We can safely use 'VB.head' here since the length of 'xs' must be larger
    -- than 'samplesMinDiagonal'.
    dimState = VS.length $ VB.head xs'
    dimMs = L.rows $ L.unSym ms
    (_, ss, xsNormalized) = S.scale xs''
    sigmaNormalized = L.unSym $ either error fst $ S.graphicalLasso 0.1 xsNormalized
    sigma = S.rescaleSWith ss sigmaNormalized
    msNew = L.inv sigma

hTuningFunctionWith ::
  Dimension ->
  -- Conversion from value to vector.
  (a -> Positions) ->
  HTuningConf ->
  Maybe (TuningFunction a)
hTuningFunctionWith n toVec (HTuningConf lc mc) = case (lc, mc) of
  (HNoTuneLeapfrog, HNoTuneMasses) -> Nothing
  (_, _) -> Just $
    \tt pdim ar xs (_, ts) ->
      let (HParamsI eps la ms tpv tpf hd) =
            -- NOTE: Use error here, because a dimension mismatch is a serious bug.
            either error id $ fromAuxiliaryTuningParameters n ts
          (TParamsVar epsMean h m) = tpv
          (TParamsFixed eps0 mu ga t0 ka) = tpf
          ms' = case mc of
            HNoTuneMasses -> ms
            HTuneDiagonalMassesOnly -> tuneDiagonalMassesOnly toVec xs ms
            HTuneAllMasses -> tuneAllMasses toVec xs ms
          hd' = case mc of
            HNoTuneMasses -> hd
            _ -> getHData ms'
          (eps'', epsMean'', h'') = case lc of
            HNoTuneLeapfrog -> (eps, epsMean, h)
            HTuneLeapfrog ->
              let delta = getOptimalRate pdim
                  c = recip $ m + t0
                  h' = (1.0 - c) * h + c * (delta - ar)
                  logEps' = mu - (sqrt m / ga) * h'
                  eps' = exp logEps'
                  mMKa = m ** (negate ka)
                  epsMean' = exp $ mMKa * logEps' + (1 - mMKa) * log epsMean
               in (eps', epsMean', h')
          eps''' = case tt of
            NormalTuningStep -> eps''
            LastTuningStep -> epsMean''
          tpv' = TParamsVar epsMean'' h'' (m + 1.0)
       in (eps''' / eps0, toAuxiliaryTuningParameters $ HParamsI eps''' la ms' tpv' tpf hd')

checkHStructureWith :: Foldable s => Masses -> HStructure s -> Maybe String
checkHStructureWith ms (HStructure x toVec fromVec)
  | toList (fromVec x xVec) /= toList x = eWith "'fromVectorWith x (toVector x) /= x' for sample state."
  | L.size xVec /= nrows = eWith "Mass matrix and 'toVector x' have different sizes for sample state."
  | otherwise = Nothing
  where
    eWith m = Just $ "checkHStructureWith: " <> m
    nrows = L.rows $ L.unSym ms
    xVec = toVec x

-- Generate momenta for a new iteration.
generateMomenta ::
  StatefulGen g m =>
  Mu ->
  Masses ->
  g ->
  m Momenta
generateMomenta mu masses gen = do
  seed <- uniformM gen
  let momenta = L.gaussianSample seed 1 mu masses
  return $ L.flatten momenta

-- TODO (high): Use a sparse matrix approach for the log density of the
-- multivariate normal, similar to McmcDate.

-- Compute exponent of kinetic energy.
exponentialKineticEnergy ::
  MassesInv ->
  Momenta ->
  Log Double
exponentialKineticEnergy msInvH xs =
  Exp $ (-0.5) * ((xs L.<# msInv) L.<.> xs)
  where
    msInv = L.unSym msInvH

-- Function calculating target value and gradient.
--
-- The function acts on the subset of the state manipulated by the proposal but
-- the value and gradient have to be calculated for the complete state. The
-- reason is that parameters untouched by the Hamiltonian proposal may affect
-- the result or the gradient.
--
-- Make sure that the value is calculated lazily because many times, only the
-- gradient is required.
type Target = Positions -> (Log Double, Positions)

-- Leapfrog integrator.
leapfrog ::
  Target ->
  MassesInv ->
  --
  LeapfrogTrajectoryLength ->
  LeapfrogScalingFactor ->
  --
  Positions ->
  Momenta ->
  -- | (New positions, new momenta, old target, new target).
  --
  -- Fail if state is not valid.
  Maybe (Positions, Momenta, Log Double, Log Double)
leapfrog tF msInv l eps q p = do
  -- The first half step of the momenta.
  (x, pHalf) <-
    let (x, pHalf) = leapfrogStepMomenta (0.5 * eps) tF q p
     in if x > 0.0
          then Just (x, pHalf)
          else Nothing
  -- L-1 full steps for positions and momenta. This gives the positions q_{L-1},
  -- and the momenta p_{L-1/2}.
  (qLM1, pLM1Half) <- go (l - 1) $ Just $ (q, pHalf)
  -- The last full step of the positions.
  let qL = leapfrogStepPositions msInv eps qLM1 pLM1Half
  -- The last half step of the momenta.
  (x', pL) <-
    let (x', pL) = leapfrogStepMomenta (0.5 * eps) tF qL pLM1Half
     in if x' > 0.0
          then Just (x', pL)
          else Nothing
  return (qL, pL, x, x')
  where
    go _ Nothing = Nothing
    go n (Just (qs, ps))
      | n <= 0 = Just (qs, ps)
      | otherwise =
          let qs' = leapfrogStepPositions msInv eps qs ps
              (x, ps') = leapfrogStepMomenta eps tF qs' p
           in if x > 0.0
                then go (n - 1) $ Just $ (qs', ps')
                else Nothing

leapfrogStepMomenta ::
  LeapfrogScalingFactor ->
  Target ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New momenta; also return value target function to be collected at the end
  -- of the leapfrog integration.
  (Log Double, Momenta)
leapfrogStepMomenta eps tf q p = (x, p + L.scale eps g)
  where
    (x, g) = tf q

leapfrogStepPositions ::
  MassesInv ->
  LeapfrogScalingFactor ->
  -- Current positions.
  Positions ->
  -- Current momenta.
  Momenta ->
  -- New positions.
  Positions
leapfrogStepPositions msInv eps q p = q + (L.unSym msInvEps L.#> p)
  where
    msInvEps = L.scale eps msInv
