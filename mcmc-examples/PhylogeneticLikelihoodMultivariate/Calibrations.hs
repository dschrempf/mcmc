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

import ELynx.Tree
import Mcmc.Tree
import Numeric.Log

-- | A calibration is specified by a name, a node at given path, and height boundaries.
--
-- For example, @let c = ("MyRootCalibration", [], YOUNG, OLD)@ ensures that the
-- root node is older than YOUNG, and younger than OLD.
type Calibration = (String, Path, Interval)

-- | Calibration prior with uniform soft bounds.
--
-- For a given set of calibrations, the absolute height of the time tree, and
-- the relative time tree, calculate the calibration prior.
--
-- The calibrations have to be precomputed with 'getCalibrations'. The reason is
-- that finding the nodes on the tree is a slow process not to be repeated after
-- each proposal.
calibrations :: HasHeight a => [Calibration] -> Double -> Tree e a -> [Log Double]
calibrations xs h t =
  [calibrateUniformSoft 1e-3 (transformInterval (recip h) i) x t | (_, x, i) <- xs]

-- | Find and calibrate the calibrated nodes on the tree.
getCalibrations :: Tree e Name -> [Calibration]
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
    cladeBrachypodium t
  ]

cladeRoot :: Tree e Name -> Calibration
cladeRoot t = ("cladeRoot", mrcaUnsafe ["Uronema_belka", "Brachypodium_distachyon"] t, properInterval 940.4 1891)

cladeChlorophyceae :: Tree e Name -> Calibration
cladeChlorophyceae t = ("cladeChlorophyceae", mrcaUnsafe ["Uronema_belka", "Monomastix_opisthostigma"] t, properInterval 940.4 1891)

cladeStreptophyta :: Tree e Name -> Calibration
cladeStreptophyta t = ("cladeStreptophyta", mrcaUnsafe ["Spirotaenia_minuta", "Brachypodium_distachyon"] t, properInterval 469 1891)

cladeEmbryophyta :: Tree e Name -> Calibration
cladeEmbryophyta t = ("cladeEmbryophyta", mrcaUnsafe ["Marchantia_polymorpha", "Brachypodium_distachyon"] t, properInterval 469 515.5)

cladeSetaphyta :: Tree e Name -> Calibration
cladeSetaphyta t = ("cladeSetaphyta", mrcaUnsafe ["Marchantia_polymorpha", "Takakia_lepidozioides"] t, properInterval 381.7 515.5)

cladePelliidae :: Tree e Name -> Calibration
cladePelliidae t = ("cladePelliidae", mrcaUnsafe ["Pellia_neesiana", "Pallavicinia_lyellii"] t, properInterval 160.4 515.5)

cladeMarchantiales :: Tree e Name -> Calibration
cladeMarchantiales t = ("cladeMarchantiales", mrcaUnsafe ["Marchantia_polymorpha", "Sphaerocarpos_texanus"] t, properInterval 227 515.5)

cladeRicciales :: Tree e Name -> Calibration
cladeRicciales t = ("cladeRicciales", mrcaUnsafe ["Riccia_berychiana", "Conocephalum_conicum"] t, properInterval 227 515.5)

cladeJungermanniidae :: Tree e Name -> Calibration
cladeJungermanniidae t = ("cladeJungermanniidae", mrcaUnsafe ["Porella_pinnata", "Metzgeria_crassipilis"] t, properInterval 144.9 515.5)

cladeJungermanniales :: Tree e Name -> Calibration
cladeJungermanniales t = ("cladeJungermanniales", mrcaUnsafe ["Porella_pinnata", "Schistochila_sp"] t, properInterval 112.7 515.5)

cladePorellineae :: Tree e Name -> Calibration
cladePorellineae t = ("cladePorellineae", mrcaUnsafe ["Porella_pinnata", "Lejeuneaceae_sp"] t, properInterval 98.17 515.5)

cladeRadulaceae :: Tree e Name -> Calibration
cladeRadulaceae t = ("cladeRadulaceae", mrcaUnsafe ["Radula_lindenbergia", "Lejeuneaceae_sp"] t, properInterval 98.17 515.5)

cladeFrullaniaceae :: Tree e Name -> Calibration
cladeFrullaniaceae t = ("cladeFrullaniaceae", mrcaUnsafe ["Frullania", "Lejeuneaceae_sp"] t, properInterval 98.17 515.5)

cladeSphagnopsida :: Tree e Name -> Calibration
cladeSphagnopsida t = ("cladeSphagnopsida", mrcaUnsafe ["Sphagnum_recurvatum", "Polytrichum_commune"] t, properInterval 330.7 515.5)

cladePolytrichopsida :: Tree e Name -> Calibration
cladePolytrichopsida t = ("cladePolytrichopsida", mrcaUnsafe ["Polytrichum_commune", "Ceratodon_purpureus"] t, properInterval 271.8 515.5)

cladePolytrichaceae :: Tree e Name -> Calibration
cladePolytrichaceae t = ("cladePolytrichaceae", mrcaUnsafe ["Polytrichum_commune", "Tetraphis_pellucida"] t, properInterval 133.3 515.5)

cladePolytrichum :: Tree e Name -> Calibration
cladePolytrichum t = ("cladePolytrichum", mrcaUnsafe ["Polytrichum_commune", "Atrichum_angustatum"] t, properInterval 82.9 515.5)

cladeFunariidae :: Tree e Name -> Calibration
cladeFunariidae t = ("cladeFunariidae", mrcaUnsafe ["Physcomitrella_patens", "Ceratodon_purpureus"] t, properInterval 268.3 515.5)

cladeDicraniidae :: Tree e Name -> Calibration
cladeDicraniidae t = ("cladeDicraniidae", mrcaUnsafe ["Thuidium_delicatulum", "Ceratodon_purpureus"] t, properInterval 133.3 515.5)

cladeHypnanae :: Tree e Name -> Calibration
cladeHypnanae t = ("cladeHypnanae", mrcaUnsafe ["Bryum_argenteum", "Brachypodium_distachyon"] t, properInterval 133.3 515.5)

cladeTracheophyta :: Tree e Name -> Calibration
cladeTracheophyta t = ("cladeTracheophyta", mrcaUnsafe ["Isoetes_sp", "Brachypodium_distachyon"] t, properInterval 420.7 451)

cladeLycopodiopphyta :: Tree e Name -> Calibration
cladeLycopodiopphyta t = ("cladeLycopodiopphyta", mrcaUnsafe ["Isoetes_sp", "Huperzia_lucidula"] t, properInterval 392.1 451)

cladeIsoetales :: Tree e Name -> Calibration
cladeIsoetales t = ("cladeIsoetales", mrcaUnsafe ["Isoetes_sp", "Selaginella_kraussiana"] t, properInterval 386.8 451)

cladeSelaginellaceae :: Tree e Name -> Calibration
cladeSelaginellaceae t = ("cladeSelaginellaceae", mrcaUnsafe ["Selaginella_selaginoides", "Selaginella_kraussiana"] t, properInterval 323.8 451)

cladeStachygynandrum :: Tree e Name -> Calibration
cladeStachygynandrum t = ("cladeStachygynandrum", mrcaUnsafe ["Selaginella_apoda", "Selaginella_kraussiana"] t, properInterval 98.17 451)

cladeLycopodioideae :: Tree e Name -> Calibration
cladeLycopodioideae t = ("cladeLycopodioideae", mrcaUnsafe ["Pseudolycopodiella_caroliniana", "Lycopodium_deuterodensum"] t, properInterval 199 451)

cladeEuphyllophyta :: Tree e Name -> Calibration
cladeEuphyllophyta t = ("cladeEuphyllophyta", mrcaUnsafe ["Equisetum_hymale", "Brachypodium_distachyon"] t, properInterval 385.5 451)

cladeMonilophyta :: Tree e Name -> Calibration
cladeMonilophyta t = ("cladeMonilophyta", mrcaUnsafe ["Equisetum_hymale", "Psilotum_nudum"] t, properInterval 384.7 451)

cladeEquisetum :: Tree e Name -> Calibration
cladeEquisetum t = ("cladeEquisetum", mrcaUnsafe ["Equisetum_hymale", "Equisetum_diffusum"] t, properInterval 64.96 451)

cladeMarattiales :: Tree e Name -> Calibration
cladeMarattiales t = ("cladeMarattiales", mrcaUnsafe ["Marattia_attenuata", "Danaea_nodosa"] t, properInterval 176 451)

cladeEusporangiates :: Tree e Name -> Calibration
cladeEusporangiates t = ("cladeEusporangiates", mrcaUnsafe ["Marattia_attenuata", "Dipteris_conjugata"] t, properInterval 318.71 451)

cladeGleicheniales :: Tree e Name -> Calibration
cladeGleicheniales t = ("cladeGleicheniales", mrcaUnsafe ["Dipteris_conjugata", "Thyrsopteris_elegans"] t, properInterval 268.3 451)

cladeCyatheales :: Tree e Name -> Calibration
cladeCyatheales t = ("cladeCyatheales", mrcaUnsafe ["Thyrsopteris_elegans", "Polystichum_acrostichoides"] t, properInterval 178 451)

cladeLindsaceae :: Tree e Name -> Calibration
cladeLindsaceae t = ("cladeLindsaceae", mrcaUnsafe ["Lindsaea_linearis", "Polystichum_acrostichoides"] t, properInterval 100.5 451)

cladeCystodiaceae :: Tree e Name -> Calibration
cladeCystodiaceae t = ("cladeCystodiaceae", mrcaUnsafe ["Cystodium_sorbifolium", "Polystichum_acrostichoides"] t, properInterval 98.17 451)

cladePteridaceae :: Tree e Name -> Calibration
cladePteridaceae t = ("cladePteridaceae", mrcaUnsafe ["Pteris_vittata", "Polystichum_acrostichoides"] t, properInterval 100.1 451)

cladeEupolypods :: Tree e Name -> Calibration
cladeEupolypods t = ("cladeEupolypods", mrcaUnsafe ["Gymnocarpium_dryopteris", "Polystichum_acrostichoides"] t, properInterval 71.5 451)

cladeSpermatophyta :: Tree e Name -> Calibration
cladeSpermatophyta t = ("cladeSpermatophyta", mrcaUnsafe ["Ginkgo_biloba", "Brachypodium_distachyon"] t, properInterval 308.14 365.6)

cladeAcrogymnospermae :: Tree e Name -> Calibration
cladeAcrogymnospermae t = ("cladeAcrogymnospermae", mrcaUnsafe ["Ginkgo_biloba", "Taxus_baccata"] t, properInterval 308.14 365.6)

cladeCycadales :: Tree e Name -> Calibration
cladeCycadales t = ("cladeCycadales", mrcaUnsafe ["Cycas_micholitzii", "Ginkgo_biloba"] t, properInterval 264.7 365.6)

cladeGnetum :: Tree e Name -> Calibration
cladeGnetum t = ("cladeGnetum", mrcaUnsafe ["Gnetum_montanum", "Ephedra_sinica"] t, properInterval 110 321.4)

cladePinopsida :: Tree e Name -> Calibration
cladePinopsida t = ("cladePinopsida", mrcaUnsafe ["Gnetum_montanum", "Cedrus_libani"] t, properInterval 153.6 321.4)

cladePinaceae :: Tree e Name -> Calibration
cladePinaceae t = ("cladePinaceae", mrcaUnsafe ["Pinus_parviflora", "Cedrus_libani"] t, properInterval 129 321.4)

cladeParviflora :: Tree e Name -> Calibration
cladeParviflora t = ("cladeParviflora", mrcaUnsafe ["Pinus_parviflora", "Pinus_radiata"] t, properInterval 89 321.4)

cladeRadiata :: Tree e Name -> Calibration
cladeRadiata t = ("cladeRadiata", mrcaUnsafe ["Pinus_jeffreyi", "Pinus_radiata"] t, properInterval 12 321.4)

cladePonderosa :: Tree e Name -> Calibration
cladePonderosa t = ("cladePonderosa", mrcaUnsafe ["Pinus_jeffreyi", "Pinus_ponderosa"] t, properInterval 6 321.4)

cladeTaxaceae :: Tree e Name -> Calibration
cladeTaxaceae t = ("cladeTaxaceae", mrcaUnsafe ["Taxus_baccata", "Juniperus_scopulorum"] t, properInterval 201 321.4)

cladeJuniperus :: Tree e Name -> Calibration
cladeJuniperus t = ("cladeJuniperus", mrcaUnsafe ["Cunninghamia_lanceolata", "Juniperus_scopulorum"] t, properInterval 83 321.4)

cladeAngiospermae :: Tree e Name -> Calibration
cladeAngiospermae t = ("cladeAngiospermae", mrcaUnsafe ["Amborella_trichopoda", "Brachypodium_distachyon"] t, properInterval 125 247.0)

cladeNympheales :: Tree e Name -> Calibration
cladeNympheales t = ("cladeNympheales", mrcaUnsafe ["Nuphar_advena", "Brachypodium_distachyon"] t, properInterval 125 247.0)

cladeAustrobaileyales :: Tree e Name -> Calibration
cladeAustrobaileyales t = ("cladeAustrobaileyales", mrcaUnsafe ["Illicium_parviflorum", "Brachypodium_distachyon"] t, properInterval 125 247.0)

cladeMesangiospermae :: Tree e Name -> Calibration
cladeMesangiospermae t = ("cladeMesangiospermae", mrcaUnsafe ["Sarcandra_glabra", "Brachypodium_distachyon"] t, properInterval 125 247.0)

cladeMagnoliids :: Tree e Name -> Calibration
cladeMagnoliids t = ("cladeMagnoliids", mrcaUnsafe ["Persea_borbonia", "Saruma_henryi"] t, properInterval 110.8 247.0)

cladePiperales :: Tree e Name -> Calibration
cladePiperales t = ("cladePiperales", mrcaUnsafe ["Houttuynia_cordata", "Saruma_henryi"] t, properInterval 44.3 247.0)

cladeEudicots :: Tree e Name -> Calibration
cladeEudicots t = ("cladeEudicots", mrcaUnsafe ["Podophyllum_peltatum", "Ipomoea_purpurea"] t, properInterval 119.6 128.63)

cladeVitales :: Tree e Name -> Calibration
cladeVitales t = ("cladeVitales", mrcaUnsafe ["Vitis_vinifera", "Ipomoea_purpurea"] t, properInterval 85.8 128.63)

cladeRosids :: Tree e Name -> Calibration
cladeRosids t = ("cladeRosids", mrcaUnsafe ["Kochia_scoparia", "Ipomoea_purpurea"] t, properInterval 85.8 128.63)

cladeEricales :: Tree e Name -> Calibration
cladeEricales t = ("cladeEricales", mrcaUnsafe ["Diospyros_malabarica", "Ipomoea_purpurea"] t, properInterval 85.8 128.63)

cladeMyrtales :: Tree e Name -> Calibration
cladeMyrtales t = ("cladeMyrtales", mrcaUnsafe ["Larrea_tridentata", "Oenothera_rosea"] t, properInterval 83.3 128.63)

cladeAsteraceae :: Tree e Name -> Calibration
cladeAsteraceae t = ("cladeAsteraceae", mrcaUnsafe ["Tanacetum_parthenium", "Catharanthus_roseus"] t, properInterval 41.5 128.63)

cladeSalicaceae :: Tree e Name -> Calibration
cladeSalicaceae t = ("cladeSalicaceae", mrcaUnsafe ["Populus_trichocarpa", "Hibiscus_cannabinus"] t, properInterval 48.57 128.63)

cladeSolanales :: Tree e Name -> Calibration
cladeSolanales t = ("cladeSolanales", mrcaUnsafe ["Solanum_tuberosum", "Ipomoea_purpurea"] t, properInterval 37.3 128.63)

cladeMonocots :: Tree e Name -> Calibration
cladeMonocots t = ("cladeMonocots", mrcaUnsafe ["Acorus_americanus", "Brachypodium_distachyon"] t, properInterval 119.5 128.63)

cladeDioscoreales :: Tree e Name -> Calibration
cladeDioscoreales t = ("cladeDioscoreales", mrcaUnsafe ["Dioscorea_villosa", "Brachypodium_distachyon"] t, properInterval 119.5 128.63)

cladeRiponogaceae :: Tree e Name -> Calibration
cladeRiponogaceae t = ("cladeRiponogaceae", mrcaUnsafe ["Smilax_bona_nox", "Colchicum_autumnale"] t, properInterval 50.5 128.63)

cladeArecales :: Tree e Name -> Calibration
cladeArecales t = ("cladeArecales", mrcaUnsafe ["Sabal_bermudana", "Brachypodium_distachyon"] t, properInterval 83.41 128.63)

cladePoales :: Tree e Name -> Calibration
cladePoales t = ("cladePoales", mrcaUnsafe ["Zea_mays", "Brachypodium_distachyon"] t, properInterval 66 128.63)

cladeBrachypodium :: Tree e Name -> Calibration
cladeBrachypodium t = ("cladeBrachypodium", mrcaUnsafe ["Oryza_sativa", "Brachypodium_distachyon"] t, properInterval 33.7 128.63)
