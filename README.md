# LocalAdaption_Nonbreeding_BGP
This repository contains code to perform analyses described in the following paper:

Migratory connectivity and ecological distinctiveness shape local adaptation to nonbreeding environments in migratory birds

Sheela P. Turbek, Christen M. Bossu, Matthew DeSaix, Marina Rodriguez, Joan Ferrer Obiol, Christine Rayne, Eric C. Anderson, Nicholas Bayly, Eben H. Paxton, Ana M. González, Darshan Narang, Thomas B. Smith, Marius Somveille, Sergio Gómez Villaverde, Mary J. Whitfield, Kevin Winker, and Kristen C. Ruegg

# Abstract

Adaptation to local environmental conditions plays a key role in the generation of biodiversity. Nonetheless, quantifying local adaptation during the nonbreeding season remains challenging in migratory species, which experience heterogeneous environments throughout the year. Theory predicts that strong migratory connectivity should promote local adaptation to the nonbreeding environment via reduced gene flow between individuals overwintering in ecologically distinct regions; however, this prediction has not been tested. Using genomic data spanning the breeding and nonbreeding ranges of three migratory songbirds, we demonstrate that stronger migratory connectivity and more distinct nonbreeding climatic niches are linked to stronger signals of local adaptation during the nonbreeding period. Additionally, we identify numerous climate-associated loci located near genes repeatedly detected across comparisons, highlighting candidate genes for climate adaptation in birds. This study provides some of the first clear evidence that seasonal migration can shape adaptive divergence during the nonbreeding season.

# Code
Investigating local adaptation on the nonbreeding ground of 3 migratory species. This includes investigating species and node-specific network migratory connectivity (NMC_XY) or the spread of the breeding populations on the nonbreeding ground, the environmental distinctiveness of the populations on the nonbreeding ground and GEA analyses using RDA.


1a.RDA_AMRE.Rmd : Code to replicate RDA analyses for American redstart breeding and nonbreeding populations. Genomic data can be downloaded from dryad.

1b.RDA_SWTH.Rmd : Code to replicate RDA analyses for Swainson's thrush breeding and nonbreeding populations. Genomic data can be downloaded from dryad.

1c.RDA_WIFL.Rmd : Code to replicate RDA analyses for Willow flycatcher breeding and nonbreeding populations. Genomic data can be downloaded from dryad.

2.Climate_distinctiveness_nonbreeding.R: Code to replicate the environmental distinctiveness analyses. Large rasters and shapefiles will need to be downloaded separately, however, urls to their location are provided.

3a.MigConnectivity_AMRE.Rmd : Code and data to replicate the MigConnectivity calculation of NMC for American redstart.

3b.MigConnectivity_SWTH.Rmd : Code and data to replicate the MigConnectivity calculation of NMC for Swainson's thrush.

3c.MigConnectivity_WIFL.Rmd : Code and data to replicate the MigConnectivity calculation of NMC for Willow flycatcher.

# Data
The genomic data including allele frequencies for breeding and nonbreeding populations can be found on dryad: XXX.
