{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Calibrations
-- Description :  Calibrations from fossil data
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Aug  3 22:37:27 2020.
module Calibrations
  ( Calibration,
    calibrations,
    getCalibrations,
  )
where

import qualified Data.ByteString.Char8 as BS
import Data.Maybe
import ELynx.Tree
import Mcmc.Tree
import Numeric.Log

-- | Calibrate a node with given path at given age.
type Calibration = (Path, Double, Double)

-- | Calibration prior with uniform soft bounds.
--
-- For a given set of calibrations, the absolute height of the time tree, and
-- the relative time tree, calculate the calibration prior.
--
-- The calibrations have to be pre-computed with 'getCalibrations'. The reason
-- is that finding the nodes on the tree is a slow process that should not be
-- repeated.
calibrations :: [Calibration] -> Double -> Tree Double Double -> [Log Double]
calibrations xs h t =
  [calibrateUniformSoft 1e-3 (a / h) (b / h) x t | (x, a, b) <- xs]

-- | Find and calibrate the calibrated nodes on the tree.
getCalibrations :: Tree e BS.ByteString -> [Calibration]
getCalibrations t =
  [ cladeRoot t,
    cladeChlorophyceae t,
    cladeStreptophyta t,
    cladeEmbryophyta t,
    cladeSetaphyta t,
    cladePelliidae t,
    cladeMarchantiales t,
    cladeRicciales t,
    cladeJungermanniidae t,
    cladeJungermanniales t,
    cladePorellineae t,
    cladeRadulaceae t,
    cladeFrullaniaceae t,
    cladeSphagnopsida t,
    cladePolytrichopsida t,
    cladePolytrichaceae t,
    cladePolytrichum t,
    cladeFunariidae t,
    cladeDicraniidae t,
    cladeHypnanae t,
    cladeTracheophyta t,
    cladeLycopodiopphyta t,
    cladeIsoetales t,
    cladeSelaginellaceae t,
    cladeStachygynandrum t,
    cladeLycopodioideae t,
    cladeEuphyllophyta t,
    cladeMonilophyta t,
    cladeEquisetum t,
    cladeMarattiales t,
    cladeEusporangiates t,
    cladeGleicheniales t,
    cladeCyatheales t,
    cladeLindsaceae t,
    cladeCystodiaceae t,
    cladePteridaceae t,
    cladeEupolypods t,
    cladeSpermatophyta t,
    cladeAcrogymnospermae t,
    cladeCycadales t,
    cladeGnetum t,
    cladePinopsida t,
    cladePinaceae t,
    cladeParviflora t,
    cladeRadiata t,
    cladePonderosa t,
    cladeTaxaceae t,
    cladeJuniperus t,
    cladeAngiospermae t,
    cladeNympheales t,
    cladeAustrobaileyales t,
    cladeMesangiospermae t,
    cladeMagnoliids t,
    cladePiperales t,
    cladeEudicots t,
    cladeVitales t,
    cladeRosids t,
    cladeEricales t,
    cladeMyrtales t,
    cladeAsteraceae t,
    cladeSalicaceae t,
    cladeSolanales t,
    cladeMonocots t,
    cladeDioscoreales t,
    cladeRiponogaceae t,
    cladeArecales t,
    cladePoales t,
    cladeBrachypodium t]

mrca' :: [BS.ByteString] -> Tree e BS.ByteString -> Path
mrca' xs = fromMaybe (error $ "mrca': Could not get MRCA for: " <> show xs) . mrca xs

cladeRoot :: Tree e BS.ByteString -> (Path, Double, Double)
cladeChlorophyceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeStreptophyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladeEmbryophyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladeSetaphyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladePelliidae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeMarchantiales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeRicciales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeJungermanniidae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeJungermanniales :: Tree e BS.ByteString -> (Path, Double, Double)
cladePorellineae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeRadulaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeFrullaniaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeSphagnopsida :: Tree e BS.ByteString -> (Path, Double, Double)
cladePolytrichopsida :: Tree e BS.ByteString -> (Path, Double, Double)
cladePolytrichaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladePolytrichum :: Tree e BS.ByteString -> (Path, Double, Double)
cladeFunariidae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeDicraniidae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeHypnanae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeTracheophyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladeLycopodiopphyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladeIsoetales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeSelaginellaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeStachygynandrum :: Tree e BS.ByteString -> (Path, Double, Double)
cladeLycopodioideae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeEuphyllophyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladeMonilophyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladeEquisetum :: Tree e BS.ByteString -> (Path, Double, Double)
cladeMarattiales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeEusporangiates :: Tree e BS.ByteString -> (Path, Double, Double)
cladeGleicheniales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeCyatheales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeLindsaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeCystodiaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladePteridaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeEupolypods :: Tree e BS.ByteString -> (Path, Double, Double)
cladeSpermatophyta :: Tree e BS.ByteString -> (Path, Double, Double)
cladeAcrogymnospermae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeCycadales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeGnetum :: Tree e BS.ByteString -> (Path, Double, Double)
cladePinopsida :: Tree e BS.ByteString -> (Path, Double, Double)
cladePinaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeParviflora :: Tree e BS.ByteString -> (Path, Double, Double)
cladeRadiata :: Tree e BS.ByteString -> (Path, Double, Double)
cladePonderosa :: Tree e BS.ByteString -> (Path, Double, Double)
cladeTaxaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeJuniperus :: Tree e BS.ByteString -> (Path, Double, Double)
cladeAngiospermae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeNympheales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeAustrobaileyales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeMesangiospermae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeMagnoliids :: Tree e BS.ByteString -> (Path, Double, Double)
cladePiperales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeEudicots :: Tree e BS.ByteString -> (Path, Double, Double)
cladeVitales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeRosids :: Tree e BS.ByteString -> (Path, Double, Double)
cladeEricales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeMyrtales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeAsteraceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeSalicaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeSolanales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeMonocots :: Tree e BS.ByteString -> (Path, Double, Double)
cladeDioscoreales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeRiponogaceae :: Tree e BS.ByteString -> (Path, Double, Double)
cladeArecales :: Tree e BS.ByteString -> (Path, Double, Double)
cladePoales :: Tree e BS.ByteString -> (Path, Double, Double)
cladeBrachypodium :: Tree e BS.ByteString -> (Path, Double, Double)

cladeRoot t = (mrca' ["Uronema_belka", "Brachypodium_distachyon"] t, 940.4, 1891)
cladeChlorophyceae t = (mrca' ["Uronema_belka", "Monomastix_opisthostigma"] t, 940.4, 1891)
cladeStreptophyta t = (mrca' ["Spirotaenia_minuta", "Brachypodium_distachyon"] t, 469, 1891)
cladeEmbryophyta t = (mrca' ["Marchantia_polymorpha", "Brachypodium_distachyon"] t, 469, 515.5)
cladeSetaphyta t = (mrca' ["Marchantia_polymorpha", "Takakia_lepidozioides"] t, 381.7, 515.5)
cladePelliidae t = (mrca' ["Pellia_neesiana", "Pallavicinia_lyellii"] t, 160.4, 515.5)
cladeMarchantiales t = (mrca' ["Marchantia_polymorpha", "Sphaerocarpos_texanus"] t, 227, 515.5)
cladeRicciales t = (mrca' ["Riccia_berychiana", "Conocephalum_conicum"] t, 227, 515.5)
cladeJungermanniidae t = (mrca' ["Porella_pinnata", "Metzgeria_crassipilis"] t, 144.9, 515.5)
cladeJungermanniales t = (mrca' ["Porella_pinnata", "Schistochila_sp"] t, 112.7, 515.5)
cladePorellineae t = (mrca' ["Porella_pinnata", "Lejeuneaceae_sp"] t, 98.17, 515.5)
cladeRadulaceae t = (mrca' ["Radula_lindenbergia", "Lejeuneaceae_sp"] t, 98.17, 515.5)
cladeFrullaniaceae t = (mrca' ["Frullania", "Lejeuneaceae_sp"] t, 98.17, 515.5)
cladeSphagnopsida t = (mrca' ["Sphagnum_recurvatum", "Polytrichum_commune"] t, 330.7, 515.5)
cladePolytrichopsida t = (mrca' ["Polytrichum_commune", "Ceratodon_purpureus"] t, 271.8, 515.5)
cladePolytrichaceae t = (mrca' ["Polytrichum_commune", "Tetraphis_pellucida"] t, 133.3, 515.5)
cladePolytrichum t = (mrca' ["Polytrichum_commune", "Atrichum_angustatum"] t, 82.9, 515.5)
cladeFunariidae t = (mrca' ["Physcomitrella_patens", "Ceratodon_purpureus"] t, 268.3, 515.5)
cladeDicraniidae t = (mrca' ["Thuidium_delicatulum", "Ceratodon_purpureus"] t, 133.3, 515.5)
cladeHypnanae t = (mrca' ["Bryum_argenteum", "Brachypodium_distachyon"] t, 133.3, 515.5)
cladeTracheophyta t = (mrca' ["Isoetes_sp", "Brachypodium_distachyon"] t, 420.7, 451)
cladeLycopodiopphyta t = (mrca' ["Isoetes_sp", "Huperzia_lucidula"] t, 392.1, 451)
cladeIsoetales t = (mrca' ["Isoetes_sp", "Selaginella_kraussiana"] t, 386.8, 451)
cladeSelaginellaceae t = (mrca' ["Selaginella_selaginoides", "Selaginella_kraussiana"] t, 323.8, 451)
cladeStachygynandrum t = (mrca' ["Selaginella_apoda", "Selaginella_kraussiana"] t, 98.17, 451)
cladeLycopodioideae t = (mrca' ["Pseudolycopodiella_caroliniana", "Lycopodium_deuterodensum"] t, 199, 451)
cladeEuphyllophyta t = (mrca' ["Equisetum_hymale", "Brachypodium_distachyon"] t, 385.5, 451)
cladeMonilophyta t = (mrca' ["Equisetum_hymale", "Psilotum_nudum"] t, 384.7, 451)
cladeEquisetum t = (mrca' ["Equisetum_hymale", "Equisetum_diffusum"] t, 64.96, 451)
cladeMarattiales t = (mrca' ["Marattia_attenuata", "Danaea_nodosa"] t, 176, 451)
cladeEusporangiates t = (mrca' ["Marattia_attenuata", "Dipteris_conjugata"] t, 318.71, 451)
cladeGleicheniales t = (mrca' ["Dipteris_conjugata", "Thyrsopteris_elegans"] t, 268.3, 451)
cladeCyatheales t = (mrca' ["Thyrsopteris_elegans", "Polystichum_acrostichoides"] t, 178, 451)
cladeLindsaceae t = (mrca' ["Lindsaea_linearis", "Polystichum_acrostichoides"] t, 100.5, 451)
cladeCystodiaceae t = (mrca' ["Cystodium_sorbifolium", "Polystichum_acrostichoides"] t, 98.17, 451)
cladePteridaceae t = (mrca' ["Pteris_vittata", "Polystichum_acrostichoides"] t, 100.1, 451)
cladeEupolypods t = (mrca' ["Gymnocarpium_dryopteris", "Polystichum_acrostichoides"] t, 71.5, 451)
cladeSpermatophyta t = (mrca' ["Ginkgo_biloba", "Brachypodium_distachyon"] t, 308.14, 365.6)
cladeAcrogymnospermae t = (mrca' ["Ginkgo_biloba", "Taxus_baccata"] t, 308.14, 365.6)
cladeCycadales t = (mrca' ["Cycas_micholitzii", "Ginkgo_biloba"] t, 264.7, 365.6)
cladeGnetum t = (mrca' ["Gnetum_montanum", "Ephedra_sinica"] t, 110, 321.4)
cladePinopsida t = (mrca' ["Gnetum_montanum", "Cedrus_libani"] t, 153.6, 321.4)
cladePinaceae t = (mrca' ["Pinus_parviflora", "Cedrus_libani"] t, 129, 321.4)
cladeParviflora t = (mrca' ["Pinus_parviflora", "Pinus_radiata"] t, 89, 321.4)
cladeRadiata t = (mrca' ["Pinus_jeffreyi", "Pinus_radiata"] t, 12, 321.4)
cladePonderosa t = (mrca' ["Pinus_jeffreyi", "Pinus_ponderosa"] t, 6, 321.4)
cladeTaxaceae t = (mrca' ["Taxus_baccata", "Juniperus_scopulorum"] t, 201, 321.4)
cladeJuniperus t = (mrca' ["Cunninghamia_lanceolata", "Juniperus_scopulorum"] t, 83, 321.4)
cladeAngiospermae t = (mrca' ["Amborella_trichopoda", "Brachypodium_distachyon"] t, 125, 247.0)
cladeNympheales t = (mrca' ["Nuphar_advena", "Brachypodium_distachyon"] t, 125, 247.0)
cladeAustrobaileyales t = (mrca' ["Illicium_parviflorum", "Brachypodium_distachyon"] t, 125, 247.0)
cladeMesangiospermae t = (mrca' ["Sarcandra_glabra", "Brachypodium_distachyon"] t, 125, 247.0)
cladeMagnoliids t = (mrca' ["Persea_borbonia", "Saruma_henryi"] t, 110.8, 247.0)
cladePiperales t = (mrca' ["Houttuynia_cordata", "Saruma_henryi"] t, 44.3, 247.0)
cladeEudicots t = (mrca' ["Podophyllum_peltatum", "Ipomoea_purpurea"] t, 119.6, 128.63)
cladeVitales t = (mrca' ["Vitis_vinifera", "Ipomoea_purpurea"] t, 85.8, 128.63)
cladeRosids t = (mrca' ["Kochia_scoparia", "Ipomoea_purpurea"] t, 85.8, 128.63)
cladeEricales t = (mrca' ["Diospyros_malabarica", "Ipomoea_purpurea"] t, 85.8, 128.63)
cladeMyrtales t = (mrca' ["Larrea_tridentata", "Oenothera_rosea"] t, 83.3, 128.63)
cladeAsteraceae t = (mrca' ["Tanacetum_parthenium", "Catharanthus_roseus"] t, 41.5, 128.63)
cladeSalicaceae t = (mrca' ["Populus_trichocarpa", "Hibiscus_cannabinus"] t, 48.57, 128.63)
cladeSolanales t = (mrca' ["Solanum_tuberosum", "Ipomoea_purpurea"] t, 37.3, 128.63)
cladeMonocots t = (mrca' ["Acorus_americanus", "Brachypodium_distachyon"] t, 119.5, 128.63)
cladeDioscoreales t = (mrca' ["Dioscorea_villosa", "Brachypodium_distachyon"] t, 119.5, 128.63)
cladeRiponogaceae t = (mrca' ["Smilax_bona_nox", "Colchicum_autumnale"] t, 50.5, 128.63)
cladeArecales t = (mrca' ["Sabal_bermudana", "Brachypodium_distachyon"] t, 83.41, 128.63)
cladePoales t = (mrca' ["Zea_mays", "Brachypodium_distachyon"] t, 66, 128.63)
cladeBrachypodium t = (mrca' ["Oryza_sativa", "Brachypodium_distachyon"] t, 33.7, 128.63)
