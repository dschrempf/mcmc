{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      :  Calibration
-- Description :  Calibrations from fossil data
-- Copyright   :  (c) Dominik Schrempf, 2020
-- License     :  GPL-3.0-or-later
--
-- Maintainer  :  dominik.schrempf@gmail.com
-- Stability   :  unstable
-- Portability :  portable
--
-- Creation date: Mon Aug  3 22:37:27 2020.
module Calibration
  ( Calibration,
    calibratedNodes,
  )
where

import qualified Data.ByteString.Char8 as BS
import Data.Maybe
import ELynx.Tree
import Mcmc.Tree

-- Calibrate a node with given path at given age.
type Calibration = (Path, Double, Double)

calibratedNodes :: Tree e BS.ByteString -> [Calibration]
calibratedNodes t =
  [ rootNode,
    cladeChlorophyceae t,
    cladeStreptophyta t,
    cladeEmbryophyta t,
    cladeSetaphyta t,
    cladeMarchantiophyta t,
    cladeMarchantiales t,
    cladeMetzgeriidae t,
    cladePorellinae t,
    cladeRadula t,
    cladeFrullania t,
    cladeBryophyta t,
    cladePolytrichopsida t,
    cladePolytrichaceae t,
    cladePolytrichum t,
    cladeFunariidae t,
    cladeDicraniidae t,
    cladeHypnanae t,
    cladeTracheophyta t,
    cladeLycopods t,
    cladeIsoetales t,
    cladeSelaginellales t,
    cladeLycopodioideae t,
    cladeEuphyllophyta t,
    cladePteridophyta t,
    cladeEquisetum t,
    cladeEusporangia t,
    cladeMarattia t,
    cladeGleicheniales t,
    cladeCyatheaceae t,
    cladePteris t,
    cladeEupolypods t,
    cladeSpermatophyta t,
    cladeAcrogymnosperms t,
    cladeCycadales t,
    cladePinopsida t,
    cladePinales t,
    cladePinaceae t,
    cladeAngiospermae t,
    cladeNympheales t,
    cladeAnita t,
    cladeMesangiospermae t,
    cladeLaurales t,
    cladePiperaceae t,
    cladeEudicots t,
    cladeProteales t,
    cladeVitales t,
    cladeRosids t,
    cladeEricales t,
    cladeAsteraceae t,
    cladeSolanales t,
    cladeSalicaceae t,
    cladeMonocots t,
    cladeDisocoreales t,
    cladeRiponogaceae t,
    cladeArecales t,
    cladePoaceae t,
    cladeBrachypodium t
  ]

rootNode :: Calibration
rootNode = (root, 940.4, 1891)

-- Convert long names to short ones.
toShort :: BS.ByteString -> BS.ByteString
toShort x = toShort' $ take 2 $ BS.split '_' x

toShort' :: [BS.ByteString] -> BS.ByteString
toShort' [x, y] =
  if l == 2
    then BS.take 7 x <> BS.singleton '_' <> y'
    else BS.take 2 x <> BS.singleton '_' <> y'
  where
    y' = BS.take 7 y
    l = BS.length y'
toShort' _ = error "toShort': List does not have two elements."

mrca' :: [BS.ByteString] -> Tree e BS.ByteString -> Path
mrca' xs = fromMaybe (error $ "mrca': Could not get MRCA for: " <> show xs) . mrca (map toShort xs)

cladeChlorophyceae :: Tree e BS.ByteString -> Calibration
cladeChlorophyceae t = (mrca' ["Uronema_belka", "Monomastix_opisthostigma"] t, 940.4, 1891)

cladeStreptophyta :: Tree e BS.ByteString -> Calibration
cladeStreptophyta t = (mrca' ["Entransia_fimbriat", "Lupinus_polyphyllus"] t, 469, 1891)

cladeEmbryophyta :: Tree e BS.ByteString -> Calibration
cladeEmbryophyta t = (mrca' ["Anthoceros_punctatus", "Lupinus_polyphyllus"] t, 469, 515.5)

cladeSetaphyta :: Tree e BS.ByteString -> Calibration
cladeSetaphyta t = (mrca' ["Sphagnum_lescurii", "Marchantia_emarginata"] t, 405, 515.5)

cladeMarchantiophyta :: Tree e BS.ByteString -> Calibration
cladeMarchantiophyta t = (mrca' ["Marchantia_emarginata", "Pellia_neesiana"] t, 405, 515.5)

cladeMarchantiales :: Tree e BS.ByteString -> Calibration
cladeMarchantiales t = (mrca' ["Sphaerocarpos_texanus", "Marchantia_emarginata"] t, 227, 515.5)

cladeMetzgeriidae :: Tree e BS.ByteString -> Calibration
cladeMetzgeriidae t = (mrca' ["Metzgeria_crassipilis", "Pallavicinia_lyellii"] t, 368.8, 515.5)

cladePorellinae :: Tree e BS.ByteString -> Calibration
cladePorellinae t = (mrca' ["Porella_pinnata", "Frullania_sp"] t, 98.17, 515.5)

cladeRadula :: Tree e BS.ByteString -> Calibration
cladeRadula t = (mrca' ["Radula_lindenbergia", "Porella_pinnata"] t, 98.17, 515.5)

cladeFrullania :: Tree e BS.ByteString -> Calibration
cladeFrullania t = (mrca' ["Frullania_sp", "Lejeuneaceae_sp"] t, 98.17, 515.5)

cladeBryophyta :: Tree e BS.ByteString -> Calibration
cladeBryophyta t = (mrca' ["Sphagnum_lescurii", "Tetraphis_pellucida"] t, 330.7, 515.5)

cladePolytrichopsida :: Tree e BS.ByteString -> Calibration
cladePolytrichopsida t = (mrca' ["Polytrichum_commune", "Buxbaumia_aphylla"] t, 271.8, 515.5)

cladePolytrichaceae :: Tree e BS.ByteString -> Calibration
cladePolytrichaceae t = (mrca' ["Polytrichum_commune", "Tetraphis_pellucida"] t, 133.3, 515.5)

cladePolytrichum :: Tree e BS.ByteString -> Calibration
cladePolytrichum t = (mrca' ["Polytrichum_commune", "Atrichum_angustatum"] t, 82.9, 515.5)

cladeFunariidae :: Tree e BS.ByteString -> Calibration
cladeFunariidae t = (mrca' ["Physcomitrella_patens", "Timmia_austriaca"] t, 268.3, 515.5)

cladeDicraniidae :: Tree e BS.ByteString -> Calibration
cladeDicraniidae t = (mrca' ["Ceratodon_purpureus", "Thuidium_delicatulum"] t, 133.3, 515.5)

cladeHypnanae :: Tree e BS.ByteString -> Calibration
cladeHypnanae t = (mrca' ["Bryum_argenteum", "Thuidium_delicatulum"] t, 133.3, 515.5)

cladeTracheophyta :: Tree e BS.ByteString -> Calibration
cladeTracheophyta t = (mrca' ["Isoetes_sp", "Lupinus_polyphyllus"] t, 420.7, 451)

cladeLycopods :: Tree e BS.ByteString -> Calibration
cladeLycopods t = (mrca' ["Isoetes_sp", "Huperzia_squarrosa"] t, 392.1, 451)

cladeIsoetales :: Tree e BS.ByteString -> Calibration
cladeIsoetales t = (mrca' ["Isoetes_sp", "Selaginella_selaginoides"] t, 386.8, 451)

cladeSelaginellales :: Tree e BS.ByteString -> Calibration
cladeSelaginellales t = (mrca' ["Selaginella_kraussiana", "Selaginella_selaginoides"] t, 323.8, 451)

cladeLycopodioideae :: Tree e BS.ByteString -> Calibration
cladeLycopodioideae t = (mrca' ["Pseudolycopodiella_caroliniana", "Lycopodium_deuterodensum"] t, 199, 451)

cladeEuphyllophyta :: Tree e BS.ByteString -> Calibration
cladeEuphyllophyta t = (mrca' ["Equisetum_diffusum", "Lupinus_polyphyllus"] t, 385.5, 451)

cladePteridophyta :: Tree e BS.ByteString -> Calibration
cladePteridophyta t = (mrca' ["Equisetum_diffusum", "Pteris_vittata"] t, 384.7, 451)

cladeEquisetum :: Tree e BS.ByteString -> Calibration
-- XXX: This calibration was erroneous. The name was corrected.
-- cladeEquisetum t = (mrca' ["Equisetum_diffusum", "Equisetum_hyemale"] t, 64.96, 451)
cladeEquisetum t = (mrca' ["Equisetum_diffusum", "Equisetum_hymale"] t, 64.96, 451)

cladeEusporangia :: Tree e BS.ByteString -> Calibration
cladeEusporangia t = (mrca' ["Danaea_nodosa", "Pteris_vittata"] t, 318.7, 451)

cladeMarattia :: Tree e BS.ByteString -> Calibration
cladeMarattia t = (mrca' ["Danaea_nodosa", "Marattia_attenuata"] t, 176, 451)

cladeGleicheniales :: Tree e BS.ByteString -> Calibration
cladeGleicheniales t = (mrca' ["Dipteris_conjugata", "Pteris_vittata"] t, 268.3, 451)

cladeCyatheaceae :: Tree e BS.ByteString -> Calibration
cladeCyatheaceae t = (mrca' ["Cyathea_spinulosa", "Thyrsopteris_elegans"] t, 178, 451)

cladePteris :: Tree e BS.ByteString -> Calibration
cladePteris t = (mrca' ["Pteris_vittata", "Cystopteris_fragilis"] t, 100.1, 451)

cladeEupolypods :: Tree e BS.ByteString -> Calibration
cladeEupolypods t = (mrca' ["Polystichum_acrostichoides", "Cystopteris_fragilis"] t, 65.5, 451)

cladeSpermatophyta :: Tree e BS.ByteString -> Calibration
cladeSpermatophyta t = (mrca' ["Ginkgo_biloba", "Lupinus_polyphyllus"] t, 308.14, 365.63)

cladeAcrogymnosperms :: Tree e BS.ByteString -> Calibration
cladeAcrogymnosperms t = (mrca' ["Ginkgo_biloba", "Cedrus_libani"] t, 308.14, 365.63)

cladeCycadales :: Tree e BS.ByteString -> Calibration
cladeCycadales t = (mrca' ["Ginkgo_biloba", "Cycas_micholitzii"] t, 264.7, 365.63)

cladePinopsida :: Tree e BS.ByteString -> Calibration
cladePinopsida t = (mrca' ["Juniperus_scopulorum", "Cedrus_libani"] t, 147, 321.48)

cladePinales :: Tree e BS.ByteString -> Calibration
cladePinales t = (mrca' ["Gnetum_montanum", "Cedrus_libani"] t, 119.6, 321.48)

cladePinaceae :: Tree e BS.ByteString -> Calibration
cladePinaceae t = (mrca' ["Pinus_ponderosa", "Cedrus_libani"] t, 99.6, 321.48)

cladeAngiospermae :: Tree e BS.ByteString -> Calibration
cladeAngiospermae t = (mrca' ["Amborella_trichopoda", "Lupinus_polyphyllus"] t, 125, 247)

cladeNympheales :: Tree e BS.ByteString -> Calibration
cladeNympheales t = (mrca' ["Nuphar_advena", "Lupinus_polyphyllus"] t, 125, 247)

cladeAnita :: Tree e BS.ByteString -> Calibration
cladeAnita t = (mrca' ["Illicium_floridanum", "Lupinus_polyphyllus"] t, 125, 247)

cladeMesangiospermae :: Tree e BS.ByteString -> Calibration
cladeMesangiospermae t = (mrca' ["Persea_borbonia", "Lupinus_polyphyllus"] t, 125, 247)

cladeLaurales :: Tree e BS.ByteString -> Calibration
cladeLaurales t = (mrca' ["Sarcandra_glabra", "Persea_borbonia"] t, 110.9, 247)

cladePiperaceae :: Tree e BS.ByteString -> Calibration
-- XXX: This calibration was erroneous. The name was corrected.
-- cladePiperaceae t = (mrca' ["Saruma_henreyi", "Houttuynia_cordata"] t, 44.3, 247)
cladePiperaceae t = (mrca' ["Saruma_henryi", "Houttuynia_cordata"] t, 44.3, 247)

cladeEudicots :: Tree e BS.ByteString -> Calibration
cladeEudicots t = (mrca' ["Escholzia_californica", "Lupinus_polyphyllus"] t, 119.6, 128.63)

cladeProteales :: Tree e BS.ByteString -> Calibration
cladeProteales t = (mrca' ["Nelumbo_nucifera", "Lupinus_polyphyllus"] t, 107.59, 128.63)

cladeVitales :: Tree e BS.ByteString -> Calibration
cladeVitales t = (mrca' ["Vitis_vinifera", "Lupinus_polyphyllus"] t, 85.8, 128.63)

cladeRosids :: Tree e BS.ByteString -> Calibration
cladeRosids t = (mrca' ["Kochia_scoparia", "Lupinus_polyphyllus"] t, 85.8, 128.63)

cladeEricales :: Tree e BS.ByteString -> Calibration
cladeEricales t = (mrca' ["Diospyros_malabarica", "Lupinus_polyphyllus"] t, 85.8, 128.63)

cladeAsteraceae :: Tree e BS.ByteString -> Calibration
cladeAsteraceae t = (mrca' ["Inula_helenium", "Solanum_tuberosum"] t, 41.5, 128.63)

cladeSolanales :: Tree e BS.ByteString -> Calibration
cladeSolanales t = (mrca' ["Ipomoea_purpurea", "Solanum_tuberosum"] t, 37.3, 128.63)

cladeSalicaceae :: Tree e BS.ByteString -> Calibration
cladeSalicaceae t = (mrca' ["Populus_euphratica", "Hibiscus_cannabinus"] t, 48.57, 128.63)

cladeMonocots :: Tree e BS.ByteString -> Calibration
cladeMonocots t = (mrca' ["Acorus_americanus", "Zea_mays"] t, 113, 128.63)

cladeDisocoreales :: Tree e BS.ByteString -> Calibration
cladeDisocoreales t = (mrca' ["Dioscorea_villosa", "Zea_mays"] t, 110.87, 128.63)

cladeRiponogaceae :: Tree e BS.ByteString -> Calibration
cladeRiponogaceae t = (mrca' ["Smilax_bona_nox", "Colchicum_autumnale"] t, 50.5, 128.63)

cladeArecales :: Tree e BS.ByteString -> Calibration
cladeArecales t = (mrca' ["Sabal_bermudana", "Zea_mays"] t, 83.41, 128.63)

cladePoaceae :: Tree e BS.ByteString -> Calibration
cladePoaceae t = (mrca' ["Oryza_sativa", "Zea_mays"] t, 65.98, 128.63)

cladeBrachypodium :: Tree e BS.ByteString -> Calibration
cladeBrachypodium t = (mrca' ["Oryza_sativa", "Brachypodium_distachyon"] t, 33.97, 128.63)
