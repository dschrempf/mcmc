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
import ELynx.Data.Tree
import NodePrior

-- Calibrate a node with given path at given age.
type Calibration = (Path, Double, Double)

calibratedNodes :: Tree e BS.ByteString -> [Calibration]
calibratedNodes _ = [rootNode]

rootNode :: Calibration
rootNode = (root, 940.4, 1891)

-- clade_chlorophyceae = clade("Uronema_belka", "Monomastix_opisthostigma")

-- # Monitor the age of node #
-- age_chlorophyceae := tmrca(timetree, clade_chlorophyceae)

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_chlorophyceae ~ dnUniform(940.4, 1891)


-- clade_streptophyta = clade("Entransia_fimbriat", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_streptophyta  := tmrca(timetree, clade_streptophyta )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_streptophyta ~ dnUniform(469, 1891)


-- clade_embryophyta = clade("Anthoceros_punctatus", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_embryophyta := tmrca(timetree, clade_embryophyta  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_embryophyta ~ dnUniform(469, 515.5)


-- clade_setaphyta = clade("Sphagnum_lescurii", "Marchantia_emarginata")

-- # Monitor the age of node #
-- age_setaphyta := tmrca(timetree, clade_setaphyta  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_setaphyta ~ dnUniform(405, 515.5)


-- clade_marchantiophyta = clade("Marchantia_emarginata", "Pellia_neesiana")

-- # Monitor the age of node #
-- age_marchantiophyta := tmrca(timetree, clade_marchantiophyta  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_marchantiophyta ~ dnUniform(405, 515.5)


-- clade_marchantiales = clade("Sphaerocarpos_texanus", "Marchantia_emarginata")

-- # Monitor the age of node #
-- age_marchantiales := tmrca(timetree, clade_marchantiales  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_marchantiales  ~ dnUniform(227, 515.5)


-- clade_metzgeriidae = clade("Metzgeria_crassipilis", "Pallavicinia_lyellii")

-- # Monitor the age of node #
-- age_metzgeriidae     := tmrca(timetree, clade_metzgeriidae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_metzgeriidae ~ dnUniform(368.8, 515.5)


-- clade_porellinae = clade("Porella_pinnata", "Frullania_sp")

-- # Monitor the age of node #
-- age_porellinae := tmrca(timetree, clade_porellinae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_porellinae ~ dnUniform(98.17, 515.5)


-- clade_radula = clade("Radula_lindenbergia", "Porella_pinnata")

-- # Monitor the age of node #
-- age_radula := tmrca(timetree, clade_radula   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_radula ~ dnUniform(98.17, 515.5)


-- clade_frullania = clade("Frullania_sp", "Lejeuneaceae_sp")

-- # Monitor the age of node #
-- age_frullania := tmrca(timetree, clade_frullania   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_frullania ~ dnUniform(98.17, 515.5)


-- clade_bryophyta = clade("Sphagnum_lescurii", "Tetraphis_pellucida")

-- # Monitor the age of node #
-- age_bryophyta := tmrca(timetree, clade_bryophyta   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_frullania ~ dnUniform(330.7, 515.5)


-- clade_polytrichopsida = clade("Polytrichum_commune", "Buxbaumia_aphylla")

-- # Monitor the age of node #
-- age_polytrichopsida := tmrca(timetree, clade_polytrichopsida   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_polytrichopsida ~ dnUniform(271.8, 515.5)


-- clade_polytrichaceae = clade("Polytrichum_commune", "Tetraphis_pellucida")

-- # Monitor the age of node #
-- age_polytrichaceae := tmrca(timetree, clade_polytrichaceae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_polytrichaceae ~ dnUniform(133.3, 515.5)


-- clade_polytrichum = clade("Polytrichum_commune", "Atrichum_angustatum")

-- # Monitor the age of node #
-- age_polytrichum := tmrca(timetree, clade_polytrichum   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_polytrichum ~ dnUniform(82.9, 515.5)


-- clade_funariidae = clade("Physcomitrella_patens", "Timmia_austriaca")

-- # Monitor the age of node #
-- age_polytrichum := tmrca(timetree, clade_funariidae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_funariidae ~ dnUniform(268.3, 515.5)


-- clade_dicraniidae = clade("Ceratodon_purpureus", "Thuidium_delicatulum")

-- # Monitor the age of node #
-- age_dicraniidae := tmrca(timetree, clade_dicraniidae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_dicraniidae ~ dnUniform(133.3, 515.5)


-- clade_hypnanae= clade("Bryum_argenteum", "Thuidium_delicatulum")

-- # Monitor the age of node #
-- age_hypnanae := tmrca(timetree, clade_hypnanae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_hypnanae ~ dnUniform(133.3, 515.5)


-- clade_tracheophyta = clade("Isoetes_sp", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_tracheophyta := tmrca(timetree, clade_tracheophyta   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_tracheophyta ~ dnUniform(420.7, 451)


-- clade_lycopods = clade("Isoetes_sp", "Huperzia_squarrosa")
-- constraints = v(clade_lycopods)

-- # Monitor the age of node #
-- age_lycopods := tmrca(timetree, clade_lycopods   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_lycopods ~ dnUniform(392.1, 451)


-- clade_isoetales = clade("Isoetes_sp", "Selaginella_selaginoides")

-- # Monitor the age of node #
-- age_isoetales := tmrca(timetree, clade_isoetales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_isoetales ~ dnUniform(386.8, 451)


-- clade_selaginellales = clade("Selaginella_kraussiana", "Selaginella_selaginoides")

-- # Monitor the age of node #
-- age_selaginellales := tmrca(timetree, clade_selaginellales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_selaginellales ~ dnUniform(323.8, 451)


-- clade_lycopodioideae = clade("Pseudolycopodiella_caroliniana", "Lycopodium_deuterodensum")

-- # Monitor the age of node #
-- age_lycopodioideae := tmrca(timetree, clade_lycopodioideae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_lycopodioideae ~ dnUniform(199, 451)


-- clade_euphyllophyta = clade("Equisetum_diffusum", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_euphyllophyta := tmrca(timetree, clade_euphyllophyta   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_euphyllophyta ~ dnUniform(385.5, 451)


-- clade_pteridophyta = clade("Equisetum_diffusum", "Pteris_vittata")

-- # Monitor the age of node #
-- age_pteridophyta := tmrca(timetree, clade_pteridophyta   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_pteridophyta ~ dnUniform(384.7, 451)


-- clade_equisetum = clade("Equisetum_diffusum", "Equisetum_hyemale")

-- # Monitor the age of node #
-- age_equisetum := tmrca(timetree, clade_equisetum   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_equisetum ~ dnUniform(64.96, 451)



-- clade_eusporangia = clade("Danaea_nodosa", "Pteris_vittata")

-- # Monitor the age of node #
-- age_eusporangia := tmrca(timetree, clade_eusporangia   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_eusporangia ~ dnUniform(318.7, 451)


-- clade_marattia = clade("Danaea_nodosa", "Marattia_attenuata")

-- # Monitor the age of node #
-- age_marattia := tmrca(timetree, clade_marattia   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_marattia ~ dnUniform(176, 451)


-- clade_gleicheniales = clade("Dipteris_conjugata", "Pteris_vittata")

-- # Monitor the age of node #
-- age_gleicheniales := tmrca(timetree, clade_gleicheniales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_gleicheniales ~ dnUniform(268.3, 451)


-- clade_cyatheaceae = clade("Cyathea_spinulosa", "Thyrsopteris_elegans")

-- # Monitor the age of node #
-- age_cyatheaceae := tmrca(timetree, clade_cyatheaceae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_cyatheaceae ~ dnUniform(178, 451)


-- clade_pteris = clade("Pteris_vittata", "Cystopteris_fragilis")

-- # Monitor the age of node #
-- age_pteris := tmrca(timetree, clade_pteris   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_pteris ~ dnUniform(100.1, 451)


-- clade_eupolypods = clade("Polystichum_acrostichoides", "Cystopteris_fragilis")
-- constraints = v(clade_eupolypods)

-- # Monitor the age of node #
-- age_eupolypods := tmrca(timetree, clade_eupolypods   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_eupolypods ~ dnUniform(65.5, 451)


-- clade_spermatophyta = clade("Ginkgo_biloba", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_spermatophyta := tmrca(timetree, clade_spermatophyta   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_spermatophyta ~ dnUniform(308.14, 365.63)


-- clade_acrogymnosperms = clade("Ginkgo_biloba", "Cedrus_libani")

-- # Monitor the age of node #
-- age_acrogymnosperms := tmrca(timetree, clade_acrogymnosperms   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_acrogymnosperms ~ dnUniform(308.14, 365.63)


-- clade_cycadales = clade("Ginkgo_biloba", "Cycas_micholitzii")

-- # Monitor the age of node #
-- age_cycadales := tmrca(timetree, clade_cycadales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_cycadales ~ dnUniform(264.7, 365.63)


-- clade_pinopsida = clade("Juniperus_scopulorum", "Cedrus_libani")

-- # Monitor the age of node #
-- age_pinopsida := tmrca(timetree, clade_pinopsida   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_pinopsida ~ dnUniform(147, 321.48)


-- clade_pinales = clade("Gnetum_montanum", "Cedrus_libani")

-- # Monitor the age of node #
-- age_pinales := tmrca(timetree, clade_pinales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_pinales ~ dnUniform(119.6, 321.48)


-- clade_pinaceae = clade("Pinus_ponderosa", "Cedrus_libani")

-- # Monitor the age of node #
-- age_pinales := tmrca(timetree, clade_pinales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_pinales ~ dnUniform(99.6, 321.48)


-- clade_angiospermae = clade("Amborella_trichopoda", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_pinales := tmrca(timetree, clade_angiospermae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_angiospermae ~ dnUniform(125, 247)


-- clade_nympheales = clade("Nuphar_advena", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_nympheales := tmrca(timetree, clade_nympheales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_nympheales ~ dnUniform(125, 247)


-- clade_anita = clade("Illicium_floridanum", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_anita := tmrca(timetree, clade_anita   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_anita ~ dnUniform(125, 247)


-- clade_mesangiospermae = clade("Persea_borbonia", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_mesangiospermae := tmrca(timetree, clade_mesangiospermae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_mesangiospermae ~ dnUniform(125, 247)


-- clade_laurales = clade("Sarcandra_glabra", "Persea_borbonia")

-- # Monitor the age of node #
-- age_laurales := tmrca(timetree, clade_laurales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_laurales ~ dnUniform(110.9, 247)

-- clade_piperaceae = clade("Saruma_henreyi", "Houttuynia_cordata")

-- # Monitor the age of node #
-- age_laurales := tmrca(timetree, clade_piperaceae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_piperaceae ~ dnUniform(44.3, 247)


-- clade_eudicots = clade("Escholzia_californica", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_eudicots := tmrca(timetree, clade_eudicots   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_eudicots ~ dnUniform(119.6, 128.63)


-- clade_proteales = clade("Nelumbo_nucifera", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_proteales := tmrca(timetree, clade_proteales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_proteales ~ dnUniform(107.59, 128.63)

-- clade_vitales = clade("Vitis_vinifera", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_vitales := tmrca(timetree, clade_vitales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_vitales ~ dnUniform(85.8, 128.63)


-- clade_rosids = clade("Kochia_scoparia", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_rosids := tmrca(timetree, clade_rosids   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_rosids ~ dnUniform(85.8, 128.63)


-- clade_ericales = clade("Diospyros_malabarica", "Lupinus_polyphyllus")

-- # Monitor the age of node #
-- age_ericales := tmrca(timetree, clade_ericales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_ericales ~ dnUniform(85.8, 128.63)


-- clade_asteraceae = clade("Inula_helenium", "Solanum_tuberosum")

-- # Monitor the age of node #
-- age_asteraceae := tmrca(timetree, clade_asteraceae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_asteraceae ~ dnUniform(41.5, 128.63)


-- clade_solanales = clade("Ipomoea_purpurea", "Solanum_tuberosum")

-- # Monitor the age of node #
-- age_solanales := tmrca(timetree, clade_solanales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_solanales ~ dnUniform(37.3, 128.63)


-- clade_salicaceae = clade("Populus_euphratica", "Hibiscus_cannabinus")

-- # Monitor the age of node #
-- age_salicaceae := tmrca(timetree, clade_salicaceae   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_salicaceae ~ dnUniform(48.57, 128.63)


-- clade_monocots = clade("Acorus_americanus", "Zea_mays")

-- # Monitor the age of node #
-- age_monocots := tmrca(timetree, clade_monocots   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_monocots ~ dnUniform(113, 128.63)


-- clade_disocoreales = clade("Dioscorea_villosa", "Zea_mays")

-- # Monitor the age of node #
-- age_disocoreales := tmrca(timetree, clade_disocoreales   )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_disocoreales ~ dnUniform(110.87, 128.63)


-- clade_riponogaceae = clade("Smilax_bona_nox", "Colchicum_autumnale")

-- # Monitor the age of node #
-- age_riponogaceae := tmrca(timetree, clade_riponogaceae  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_riponogaceae ~ dnUniform(50.5, 128.63)


-- clade_arecales = clade("Sabal_bermudana", "Zea_mays")

-- # Monitor the age of node #
-- age_arecales := tmrca(timetree, clade_arecales  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_arecales ~ dnUniform(83.41, 128.63)


-- clade_poaceae = clade("Oryza_sativa", "Zea_mays")

-- # Monitor the age of node #
-- age_poaceae := tmrca(timetree, clade_poaceae  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_poaceae ~ dnUniform(65.98, 128.63)


-- clade_brachypodium = clade("Oryza_sativa", "Brachypodium_distachyon")

-- # Monitor the age of node #
-- age_brachypodium := tmrca(timetree, clade_brachypodium  )

-- # Specify an exponetial distribution on the age of this node #
-- obs_age_brachypodium ~ dnUniform(33.97, 128.63)
