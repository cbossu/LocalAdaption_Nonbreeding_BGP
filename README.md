# LocalAdaption_Nonbreeding_BGP
This repository contains code to perform analyses described in the following paper:

XXXX

# Abstract

Adaptation to local environmental conditions plays a key role in the generation of biodiversity. Nonetheless, quantifying local adaptation in migratory species, which experience heterogeneous environments throughout the year, remains challenging. Theory suggests that strong connections between breeding and nonbreeding populations should promote local adaptation to the nonbreeding environment; however, no study has tested this prediction. We leverage genomic data from three migratory songbirds to document a positive association between the strength of migratory connections and extent of local adaptation to climate on the nonbreeding grounds. Many climate-linked loci were located near genes detected in multiple comparisons, revealing candidate genes for climate adaptation in birds. This is one of the first studies to clearly demonstrate how seasonal migration can shape adaptive divergence during the nonbreeding season.


# Code
Investigating local adaptation on the nonbreeding ground of 3 migratory species. This includes investigating species and node-specific network migratory connectivity (NMC_XY) or the spread of the breeding populations on the nonbreeding ground, the environmental distinctiveness of the populations on the nonbreeding ground and GEA analyses using RDA.


1.RDA_of_pop_allele_freq.R : Code to replicate RDA analyses run for 3 species on breeding and nonbreeding ground using population allele frequencies estimated using a snakemake pipeline, loco-pipe.

2.Climate_distinctiveness_nonbreeding.R: Code to replicate the environmental distinctiveness analysis. Large rasters and shapefiles will need to be downloaded separately, however, url to their location are provide.

3a.MigConnectivity_AMRE.Rmd : Code and data to replicate the MigConnectivity calculation of NMC for American redstart.

3b.MigConnectivity_SWTH.Rmd : Code and data to replicate the MigConnectivity calculation of NMC for Swainson's thrush.

3c.MigConnectivity_WIFL.Rmd : Code and data to replicate the MigConnectivity calculation of NMC for Willow flycatcher.

# Data
The genotypic data including genotype probabilities for breeding and nonbreeding individuals and allele frequncies for breeding and nonbreeding populations can be found on dryad:
