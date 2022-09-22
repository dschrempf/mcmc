{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      :  Mcmc.Internal.Gamma
-- Description :  Generalized gamma function for automatic differentiation
-- Copyright   :  2021 Dominik Schrempf
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  experimental
-- Portability :  portable
--
-- Creation date: Tue Jul 13 12:53:09 2021.
--
-- The code is taken from "Numeric.SpecFunctions".
module Mcmc.Internal.SpecFunctions
  ( logGammaG,
    logFactorialG,
  )
where

import Data.Typeable
import qualified Data.Vector as VB
import Numeric.Polynomial
import Numeric.SpecFunctions
import Unsafe.Coerce

mSqrtEpsG :: RealFloat a => a
mSqrtEpsG = 1.4901161193847656e-8

mEulerMascheroniG :: RealFloat a => a
mEulerMascheroniG = 0.5772156649015328606065121

-- | Generalized version of the log gamma distribution. See
-- 'Numeric.SpecFunctions.logGamma'.
logGammaG :: (Typeable a, RealFloat a) => a -> a
logGammaG z
  | typeOf z == typeRep (Proxy :: Proxy Double) = unsafeCoerce logGamma z
  | otherwise = logGammaNonDouble z
{-# SPECIALIZE logGammaG :: Double -> Double #-}

-- See 'Numeric.SpecFunctions.logGamma'.
logGammaNonDouble :: RealFloat a => a -> a
logGammaNonDouble z
  | z <= 0 = 1 / 0
  | z < mSqrtEpsG = log (1 / z - mEulerMascheroniG)
  | z < 0.5 = lgamma1_15G z (z - 1) - log z
  | z < 1 = lgamma15_2G z (z - 1) - log z
  | z <= 1.5 = lgamma1_15G (z - 1) (z - 2)
  | z < 2 = lgamma15_2G (z - 1) (z - 2)
  | z < 15 = lgammaSmallG z
  | otherwise = lanczosApproxG z

lgamma1_15G :: RealFloat a => a -> a -> a
lgamma1_15G zm1 zm2 =
  r * y
    + r
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
  r * y
    + r
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
  r * y
    + r
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

-- | Generalized version of the log factorial function. See
-- 'Numeric.SpecFunctions.logFactorial'.
logFactorialG :: forall a b. (Integral a, RealFloat b, Typeable b) => a -> b
logFactorialG n
  | typeRep (Proxy :: Proxy b) == typeRep (Proxy :: Proxy Double) = unsafeCoerce $ logFactorial n
  | otherwise = logFactorialNonDouble n
{-# SPECIALIZE logFactorialG :: Int -> Double #-}

logFactorialNonDouble :: (Integral a, RealFloat b) => a -> b
logFactorialNonDouble n
  | n < 0 = error "logFactorialNonDouble: Negative input."
  | n <= 170 = log $ VB.unsafeIndex factorialTable (fromIntegral n)
  | n < 1500 = stirling + rx * ((1 / 12) - (1 / 360) * rx * rx)
  | otherwise = stirling + (1 / 12) * rx
  where
    stirling = (x - 0.5) * log x - x + mLnSqrt2Pi
    x = fromIntegral n + 1
    rx = recip x
{-# SPECIALIZE logFactorialNonDouble :: RealFloat a => Int -> a #-}

mLnSqrt2Pi :: RealFloat a => a
mLnSqrt2Pi = 0.9189385332046727417803297364056176398613974736377834128171
{-# INLINE mLnSqrt2Pi #-}

factorialTable :: RealFloat a => VB.Vector a
{-# NOINLINE factorialTable #-}
factorialTable =
  VB.fromListN
    171
    [ 1.0,
      1.0,
      2.0,
      6.0,
      24.0,
      120.0,
      720.0,
      5040.0,
      40320.0,
      362880.0,
      3628800.0,
      3.99168e7,
      4.790016e8,
      6.2270208e9,
      8.71782912e10,
      1.307674368e12,
      2.0922789888e13,
      3.55687428096e14,
      6.402373705728e15,
      1.21645100408832e17,
      2.43290200817664e18,
      5.109094217170944e19,
      1.1240007277776077e21,
      2.5852016738884974e22,
      6.204484017332394e23,
      1.5511210043330984e25,
      4.032914611266056e26,
      1.0888869450418352e28,
      3.0488834461171384e29,
      8.841761993739702e30,
      2.6525285981219103e32,
      8.222838654177922e33,
      2.631308369336935e35,
      8.683317618811886e36,
      2.9523279903960412e38,
      1.0333147966386144e40,
      3.719933267899012e41,
      1.3763753091226343e43,
      5.23022617466601e44,
      2.0397882081197442e46,
      8.159152832478977e47,
      3.3452526613163803e49,
      1.4050061177528798e51,
      6.041526306337383e52,
      2.6582715747884485e54,
      1.1962222086548019e56,
      5.5026221598120885e57,
      2.5862324151116818e59,
      1.2413915592536073e61,
      6.082818640342675e62,
      3.0414093201713376e64,
      1.5511187532873822e66,
      8.065817517094388e67,
      4.2748832840600255e69,
      2.308436973392414e71,
      1.2696403353658275e73,
      7.109985878048634e74,
      4.0526919504877214e76,
      2.3505613312828785e78,
      1.386831185456898e80,
      8.32098711274139e81,
      5.075802138772247e83,
      3.146997326038793e85,
      1.9826083154044399e87,
      1.2688693218588415e89,
      8.24765059208247e90,
      5.44344939077443e92,
      3.647111091818868e94,
      2.4800355424368305e96,
      1.711224524281413e98,
      1.197857166996989e100,
      8.504785885678623e101,
      6.1234458376886085e103,
      4.470115461512684e105,
      3.307885441519386e107,
      2.4809140811395396e109,
      1.88549470166605e111,
      1.4518309202828586e113,
      1.1324281178206297e115,
      8.946182130782974e116,
      7.15694570462638e118,
      5.797126020747368e120,
      4.753643337012841e122,
      3.9455239697206583e124,
      3.314240134565353e126,
      2.81710411438055e128,
      2.422709538367273e130,
      2.1077572983795275e132,
      1.8548264225739844e134,
      1.650795516090846e136,
      1.4857159644817613e138,
      1.352001527678403e140,
      1.2438414054641305e142,
      1.1567725070816416e144,
      1.087366156656743e146,
      1.0329978488239058e148,
      9.916779348709496e149,
      9.619275968248211e151,
      9.426890448883246e153,
      9.332621544394413e155,
      9.332621544394415e157,
      9.425947759838358e159,
      9.614466715035125e161,
      9.902900716486179e163,
      1.0299016745145626e166,
      1.0813967582402908e168,
      1.1462805637347082e170,
      1.2265202031961378e172,
      1.3246418194518288e174,
      1.4438595832024934e176,
      1.5882455415227428e178,
      1.7629525510902446e180,
      1.974506857221074e182,
      2.2311927486598134e184,
      2.543559733472187e186,
      2.9250936934930154e188,
      3.393108684451898e190,
      3.9699371608087206e192,
      4.68452584975429e194,
      5.574585761207606e196,
      6.689502913449126e198,
      8.094298525273443e200,
      9.875044200833601e202,
      1.214630436702533e205,
      1.5061417415111406e207,
      1.8826771768889257e209,
      2.372173242880047e211,
      3.0126600184576594e213,
      3.856204823625804e215,
      4.974504222477286e217,
      6.466855489220473e219,
      8.471580690878819e221,
      1.1182486511960041e224,
      1.4872707060906857e226,
      1.9929427461615188e228,
      2.6904727073180504e230,
      3.6590428819525483e232,
      5.012888748274991e234,
      6.917786472619488e236,
      9.615723196941088e238,
      1.3462012475717523e241,
      1.898143759076171e243,
      2.6953641378881624e245,
      3.8543707171800725e247,
      5.5502938327393044e249,
      8.047926057471992e251,
      1.1749972043909107e254,
      1.7272458904546386e256,
      2.5563239178728654e258,
      3.808922637630569e260,
      5.713383956445854e262,
      8.62720977423324e264,
      1.3113358856834524e267,
      2.0063439050956823e269,
      3.0897696138473508e271,
      4.789142901463393e273,
      7.471062926282894e275,
      1.1729568794264143e278,
      1.8532718694937346e280,
      2.946702272495038e282,
      4.714723635992061e284,
      7.590705053947218e286,
      1.2296942187394494e289,
      2.0044015765453023e291,
      3.287218585534296e293,
      5.423910666131589e295,
      9.003691705778436e297,
      1.5036165148649988e300,
      2.526075744973198e302,
      4.269068009004705e304,
      7.257415615307998e306
    ]
