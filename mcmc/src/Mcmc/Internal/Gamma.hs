-- |
-- Module      :  Mcmc.Internal.Gamma
-- Description :  Generalized gamma function for automatic differentiation
-- Copyright   :  (c) 2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Tue Jul 13 12:53:09 2021.
--
-- The code is taken from "Numeric.SpecFunctions".
module Mcmc.Internal.Gamma
  ( logGammaG,
  )
where

import Data.Typeable
import qualified Data.Vector as VB
import Numeric.Polynomial
import Numeric.SpecFunctions
import Unsafe.Coerce

mSqrtEps :: RealFloat a => a
mSqrtEps = 1.4901161193847656e-8

mEulerMascheroni :: RealFloat a => a
mEulerMascheroni = 0.5772156649015328606065121

logGammaG :: (Typeable a, RealFloat a) => a -> a
logGammaG z
  | typeOf z == typeOf (0 :: Double) = unsafeCoerce $ logGamma $ unsafeCoerce z
  | otherwise = logGammaNonDouble z
{-# SPECIALIZE logGammaG :: Double -> Double #-}

-- | See 'Numeric.SpecFunctions.logGamma'.
logGammaNonDouble :: RealFloat a => a -> a
logGammaNonDouble z
  | z <= 0 = 1 / 0
  | z < mSqrtEps = log (1 / z - mEulerMascheroni)
  | z < 0.5 = lgamma1_15G z (z - 1) - log z
  | z < 1 = lgamma15_2G z (z - 1) - log z
  | z <= 1.5 = lgamma1_15G (z - 1) (z - 2)
  | z < 2 = lgamma15_2G (z - 1) (z - 2)
  | z < 15 = lgammaSmallG z
  | otherwise = lanczosApproxG z

lgamma1_15G :: RealFloat a => a -> a -> a
lgamma1_15G zm1 zm2 =
  r * y + r
    * ( evaluatePolynomial zm1 tableLogGamma_1_15PG
          / evaluatePolynomial zm1 tableLogGamma_1_15QG
      )
  where
    r = zm1 * zm2
    y = 0.52815341949462890625

tableLogGamma_1_15PG :: RealFloat a => VB.Vector a
tableLogGamma_1_15PG =
  VB.fromList
    [ 0.490622454069039543534e-1,
      -0.969117530159521214579e-1,
      -0.414983358359495381969e0,
      -0.406567124211938417342e0,
      -0.158413586390692192217e0,
      -0.240149820648571559892e-1,
      -0.100346687696279557415e-2
    ]
{-# NOINLINE tableLogGamma_1_15PG #-}

tableLogGamma_1_15QG :: RealFloat a => VB.Vector a
tableLogGamma_1_15QG =
  VB.fromList
    [ 1,
      0.302349829846463038743e1,
      0.348739585360723852576e1,
      0.191415588274426679201e1,
      0.507137738614363510846e0,
      0.577039722690451849648e-1,
      0.195768102601107189171e-2
    ]
{-# NOINLINE tableLogGamma_1_15QG #-}

lgamma15_2G :: RealFloat a => a -> a -> a
lgamma15_2G zm1 zm2 =
  r * y + r
    * ( evaluatePolynomial (-zm2) tableLogGamma_15_2PG
          / evaluatePolynomial (-zm2) tableLogGamma_15_2QG
      )
  where
    r = zm1 * zm2
    y = 0.452017307281494140625

tableLogGamma_15_2PG :: RealFloat a => VB.Vector a
tableLogGamma_15_2PG =
  VB.fromList
    [ -0.292329721830270012337e-1,
      0.144216267757192309184e0,
      -0.142440390738631274135e0,
      0.542809694055053558157e-1,
      -0.850535976868336437746e-2,
      0.431171342679297331241e-3
    ]
{-# NOINLINE tableLogGamma_15_2PG #-}

tableLogGamma_15_2QG :: RealFloat a => VB.Vector a
tableLogGamma_15_2QG =
  VB.fromList
    [ 1,
      -0.150169356054485044494e1,
      0.846973248876495016101e0,
      -0.220095151814995745555e0,
      0.25582797155975869989e-1,
      -0.100666795539143372762e-2,
      -0.827193521891290553639e-6
    ]
{-# NOINLINE tableLogGamma_15_2QG #-}

lgammaSmallG :: RealFloat a => a -> a
lgammaSmallG = go 0
  where
    go acc z
      | z < 3 = acc + lgamma2_3G z
      | otherwise = go (acc + log zm1) zm1
      where
        zm1 = z - 1

lgamma2_3G :: RealFloat a => a -> a
lgamma2_3G z =
  r * y + r
    * ( evaluatePolynomial zm2 tableLogGamma_2_3PG
          / evaluatePolynomial zm2 tableLogGamma_2_3QG
      )
  where
    r = zm2 * (z + 1)
    zm2 = z - 2
    y = 0.158963680267333984375e0

tableLogGamma_2_3PG :: RealFloat a => VB.Vector a
tableLogGamma_2_3PG =
  VB.fromList
    [ -0.180355685678449379109e-1,
      0.25126649619989678683e-1,
      0.494103151567532234274e-1,
      0.172491608709613993966e-1,
      -0.259453563205438108893e-3,
      -0.541009869215204396339e-3,
      -0.324588649825948492091e-4
    ]
{-# NOINLINE tableLogGamma_2_3PG #-}

tableLogGamma_2_3QG :: RealFloat a => VB.Vector a
tableLogGamma_2_3QG =
  VB.fromList
    [ 1,
      0.196202987197795200688e1,
      0.148019669424231326694e1,
      0.541391432071720958364e0,
      0.988504251128010129477e-1,
      0.82130967464889339326e-2,
      0.224936291922115757597e-3,
      -0.223352763208617092964e-6
    ]
{-# NOINLINE tableLogGamma_2_3QG #-}

lanczosApproxG :: RealFloat a => a -> a
lanczosApproxG z =
  (log (z + g - 0.5) - 1) * (z - 0.5)
    + log (evalRatioG tableLanczosG z)
  where
    g = 6.024680040776729583740234375

tableLanczosG :: RealFloat a => VB.Vector (a, a)
tableLanczosG =
  VB.fromList
    [ (56906521.91347156388090791033559122686859, 0),
      (103794043.1163445451906271053616070238554, 39916800),
      (86363131.28813859145546927288977868422342, 120543840),
      (43338889.32467613834773723740590533316085, 150917976),
      (14605578.08768506808414169982791359218571, 105258076),
      (3481712.15498064590882071018964774556468, 45995730),
      (601859.6171681098786670226533699352302507, 13339535),
      (75999.29304014542649875303443598909137092, 2637558),
      (6955.999602515376140356310115515198987526, 357423),
      (449.9445569063168119446858607650988409623, 32670),
      (19.51992788247617482847860966235652136208, 1925),
      (0.5098416655656676188125178644804694509993, 66),
      (0.006061842346248906525783753964555936883222, 1)
    ]
{-# NOINLINE tableLanczosG #-}

data LG a = LG !a !a

evalRatioG :: RealFloat a => VB.Vector (a, a) -> a -> a
evalRatioG coef x
  | x > 1 = fini $ VB.foldl' stepL (LG 0 0) coef
  | otherwise = fini $ VB.foldr' stepR (LG 0 0) coef
  where
    fini (LG num den) = num / den
    stepR (a, b) (LG num den) = LG (num * x + a) (den * x + b)
    stepL (LG num den) (a, b) = LG (num * rx + a) (den * rx + b)
    rx = recip x
