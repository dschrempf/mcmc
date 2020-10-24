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
type Calibration = (String, Path, Double, Double)

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
  [calibrateUniformSoft 1e-3 (a / h) (b / h) x t | (_, x, a, b) <- xs]

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

cladeRoot :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeChlorophyceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeStreptophyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeEmbryophyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeSetaphyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePelliidae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeMarchantiales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeRicciales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeJungermanniidae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeJungermanniales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePorellineae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeRadulaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeFrullaniaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeSphagnopsida :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePolytrichopsida :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePolytrichaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePolytrichum :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeFunariidae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeDicraniidae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeHypnanae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeTracheophyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeLycopodiopphyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeIsoetales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeSelaginellaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeStachygynandrum :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeLycopodioideae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeEuphyllophyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeMonilophyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeEquisetum :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeMarattiales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeEusporangiates :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeGleicheniales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeCyatheales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeLindsaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeCystodiaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePteridaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeEupolypods :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeSpermatophyta :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeAcrogymnospermae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeCycadales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeGnetum :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePinopsida :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePinaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeParviflora :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeRadiata :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePonderosa :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeTaxaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeJuniperus :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeAngiospermae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeNympheales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeAustrobaileyales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeMesangiospermae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeMagnoliids :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePiperales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeEudicots :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeVitales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeRosids :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeEricales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeMyrtales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeAsteraceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeSalicaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeSolanales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeMonocots :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeDioscoreales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeRiponogaceae :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeArecales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladePoales :: Tree e BS.ByteString -> (String, Path, Double, Double)
cladeBrachypodium :: Tree e BS.ByteString -> (String, Path, Double, Double)

cladeRoot t = ("cladeRoot", mrca' ["Uronema_belka", "Brachypodium_distachyon"] t, 940.4, 1891)
cladeChlorophyceae t = ("cladeChlorophyceae", mrca' ["Uronema_belka", "Monomastix_opisthostigma"] t, 940.4, 1891)
cladeStreptophyta t = ("cladeStreptophyta", mrca' ["Spirotaenia_minuta", "Brachypodium_distachyon"] t, 469, 1891)
cladeEmbryophyta t = ("cladeEmbryophyta", mrca' ["Marchantia_polymorpha", "Brachypodium_distachyon"] t, 469, 515.5)
cladeSetaphyta t = ("cladeSetaphyta", mrca' ["Marchantia_polymorpha", "Takakia_lepidozioides"] t, 381.7, 515.5)
cladePelliidae t = ("cladePelliidae", mrca' ["Pellia_neesiana", "Pallavicinia_lyellii"] t, 160.4, 515.5)
cladeMarchantiales t = ("cladeMarchantiales", mrca' ["Marchantia_polymorpha", "Sphaerocarpos_texanus"] t, 227, 515.5)
cladeRicciales t = ("cladeRicciales", mrca' ["Riccia_berychiana", "Conocephalum_conicum"] t, 227, 515.5)
cladeJungermanniidae t = ("cladeJungermanniidae", mrca' ["Porella_pinnata", "Metzgeria_crassipilis"] t, 144.9, 515.5)
cladeJungermanniales t = ("cladeJungermanniales", mrca' ["Porella_pinnata", "Schistochila_sp"] t, 112.7, 515.5)
cladePorellineae t = ("cladePorellineae", mrca' ["Porella_pinnata", "Lejeuneaceae_sp"] t, 98.17, 515.5)
cladeRadulaceae t = ("cladeRadulaceae", mrca' ["Radula_lindenbergia", "Lejeuneaceae_sp"] t, 98.17, 515.5)
cladeFrullaniaceae t = ("cladeFrullaniaceae", mrca' ["Frullania", "Lejeuneaceae_sp"] t, 98.17, 515.5)
cladeSphagnopsida t = ("cladeSphagnopsida", mrca' ["Sphagnum_recurvatum", "Polytrichum_commune"] t, 330.7, 515.5)
cladePolytrichopsida t = ("cladePolytrichopsida", mrca' ["Polytrichum_commune", "Ceratodon_purpureus"] t, 271.8, 515.5)
cladePolytrichaceae t = ("cladePolytrichaceae", mrca' ["Polytrichum_commune", "Tetraphis_pellucida"] t, 133.3, 515.5)
cladePolytrichum t = ("cladePolytrichum", mrca' ["Polytrichum_commune", "Atrichum_angustatum"] t, 82.9, 515.5)
cladeFunariidae t = ("cladeFunariidae", mrca' ["Physcomitrella_patens", "Ceratodon_purpureus"] t, 268.3, 515.5)
cladeDicraniidae t = ("cladeDicraniidae", mrca' ["Thuidium_delicatulum", "Ceratodon_purpureus"] t, 133.3, 515.5)
cladeHypnanae t = ("cladeHypnanae", mrca' ["Bryum_argenteum", "Brachypodium_distachyon"] t, 133.3, 515.5)
cladeTracheophyta t = ("cladeTracheophyta", mrca' ["Isoetes_sp", "Brachypodium_distachyon"] t, 420.7, 451)
cladeLycopodiopphyta t = ("cladeLycopodiopphyta", mrca' ["Isoetes_sp", "Huperzia_lucidula"] t, 392.1, 451)
cladeIsoetales t = ("cladeIsoetales", mrca' ["Isoetes_sp", "Selaginella_kraussiana"] t, 386.8, 451)
cladeSelaginellaceae t = ("cladeSelaginellaceae", mrca' ["Selaginella_selaginoides", "Selaginella_kraussiana"] t, 323.8, 451)
cladeStachygynandrum t = ("cladeStachygynandrum", mrca' ["Selaginella_apoda", "Selaginella_kraussiana"] t, 98.17, 451)
cladeLycopodioideae t = ("cladeLycopodioideae", mrca' ["Pseudolycopodiella_caroliniana", "Lycopodium_deuterodensum"] t, 199, 451)
cladeEuphyllophyta t = ("cladeEuphyllophyta", mrca' ["Equisetum_hymale", "Brachypodium_distachyon"] t, 385.5, 451)
cladeMonilophyta t = ("cladeMonilophyta", mrca' ["Equisetum_hymale", "Psilotum_nudum"] t, 384.7, 451)
cladeEquisetum t = ("cladeEquisetum", mrca' ["Equisetum_hymale", "Equisetum_diffusum"] t, 64.96, 451)
cladeMarattiales t = ("cladeMarattiales", mrca' ["Marattia_attenuata", "Danaea_nodosa"] t, 176, 451)
cladeEusporangiates t = ("cladeEusporangiates", mrca' ["Marattia_attenuata", "Dipteris_conjugata"] t, 318.71, 451)
cladeGleicheniales t = ("cladeGleicheniales", mrca' ["Dipteris_conjugata", "Thyrsopteris_elegans"] t, 268.3, 451)
cladeCyatheales t = ("cladeCyatheales", mrca' ["Thyrsopteris_elegans", "Polystichum_acrostichoides"] t, 178, 451)
cladeLindsaceae t = ("cladeLindsaceae", mrca' ["Lindsaea_linearis", "Polystichum_acrostichoides"] t, 100.5, 451)
cladeCystodiaceae t = ("cladeCystodiaceae", mrca' ["Cystodium_sorbifolium", "Polystichum_acrostichoides"] t, 98.17, 451)
cladePteridaceae t = ("cladePteridaceae", mrca' ["Pteris_vittata", "Polystichum_acrostichoides"] t, 100.1, 451)
cladeEupolypods t = ("cladeEupolypods", mrca' ["Gymnocarpium_dryopteris", "Polystichum_acrostichoides"] t, 71.5, 451)
cladeSpermatophyta t = ("cladeSpermatophyta", mrca' ["Ginkgo_biloba", "Brachypodium_distachyon"] t, 308.14, 365.6)
cladeAcrogymnospermae t = ("cladeAcrogymnospermae", mrca' ["Ginkgo_biloba", "Taxus_baccata"] t, 308.14, 365.6)
cladeCycadales t = ("cladeCycadales", mrca' ["Cycas_micholitzii", "Ginkgo_biloba"] t, 264.7, 365.6)
cladeGnetum t = ("cladeGnetum", mrca' ["Gnetum_montanum", "Ephedra_sinica"] t, 110, 321.4)
cladePinopsida t = ("cladePinopsida", mrca' ["Gnetum_montanum", "Cedrus_libani"] t, 153.6, 321.4)
cladePinaceae t = ("cladePinaceae", mrca' ["Pinus_parviflora", "Cedrus_libani"] t, 129, 321.4)
cladeParviflora t = ("cladeParviflora", mrca' ["Pinus_parviflora", "Pinus_radiata"] t, 89, 321.4)
cladeRadiata t = ("cladeRadiata", mrca' ["Pinus_jeffreyi", "Pinus_radiata"] t, 12, 321.4)
cladePonderosa t = ("cladePonderosa", mrca' ["Pinus_jeffreyi", "Pinus_ponderosa"] t, 6, 321.4)
cladeTaxaceae t = ("cladeTaxaceae", mrca' ["Taxus_baccata", "Juniperus_scopulorum"] t, 201, 321.4)
cladeJuniperus t = ("cladeJuniperus", mrca' ["Cunninghamia_lanceolata", "Juniperus_scopulorum"] t, 83, 321.4)
cladeAngiospermae t = ("cladeAngiospermae", mrca' ["Amborella_trichopoda", "Brachypodium_distachyon"] t, 125, 247.0)
cladeNympheales t = ("cladeNympheales", mrca' ["Nuphar_advena", "Brachypodium_distachyon"] t, 125, 247.0)
cladeAustrobaileyales t = ("cladeAustrobaileyales", mrca' ["Illicium_parviflorum", "Brachypodium_distachyon"] t, 125, 247.0)
cladeMesangiospermae t = ("cladeMesangiospermae", mrca' ["Sarcandra_glabra", "Brachypodium_distachyon"] t, 125, 247.0)
cladeMagnoliids t = ("cladeMagnoliids", mrca' ["Persea_borbonia", "Saruma_henryi"] t, 110.8, 247.0)
cladePiperales t = ("cladePiperales", mrca' ["Houttuynia_cordata", "Saruma_henryi"] t, 44.3, 247.0)
cladeEudicots t = ("cladeEudicots", mrca' ["Podophyllum_peltatum", "Ipomoea_purpurea"] t, 119.6, 128.63)
cladeVitales t = ("cladeVitales", mrca' ["Vitis_vinifera", "Ipomoea_purpurea"] t, 85.8, 128.63)
cladeRosids t = ("cladeRosids", mrca' ["Kochia_scoparia", "Ipomoea_purpurea"] t, 85.8, 128.63)
cladeEricales t = ("cladeEricales", mrca' ["Diospyros_malabarica", "Ipomoea_purpurea"] t, 85.8, 128.63)
cladeMyrtales t = ("cladeMyrtales", mrca' ["Larrea_tridentata", "Oenothera_rosea"] t, 83.3, 128.63)
cladeAsteraceae t = ("cladeAsteraceae", mrca' ["Tanacetum_parthenium", "Catharanthus_roseus"] t, 41.5, 128.63)
cladeSalicaceae t = ("cladeSalicaceae", mrca' ["Populus_trichocarpa", "Hibiscus_cannabinus"] t, 48.57, 128.63)
cladeSolanales t = ("cladeSolanales", mrca' ["Solanum_tuberosum", "Ipomoea_purpurea"] t, 37.3, 128.63)
cladeMonocots t = ("cladeMonocots", mrca' ["Acorus_americanus", "Brachypodium_distachyon"] t, 119.5, 128.63)
cladeDioscoreales t = ("cladeDioscoreales", mrca' ["Dioscorea_villosa", "Brachypodium_distachyon"] t, 119.5, 128.63)
cladeRiponogaceae t = ("cladeRiponogaceae", mrca' ["Smilax_bona_nox", "Colchicum_autumnale"] t, 50.5, 128.63)
cladeArecales t = ("cladeArecales", mrca' ["Sabal_bermudana", "Brachypodium_distachyon"] t, 83.41, 128.63)
cladePoales t = ("cladePoales", mrca' ["Zea_mays", "Brachypodium_distachyon"] t, 66, 128.63)
cladeBrachypodium t = ("cladeBrachypodium", mrca' ["Oryza_sativa", "Brachypodium_distachyon"] t, 33.7, 128.63)
