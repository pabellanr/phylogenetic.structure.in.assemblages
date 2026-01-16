# phylogenetic.structure.in.assemblages
Code to assess phylogenetic clustering in an assemblage in relation to a species pool and identify which clades significantly contribute to the phylogenetic structure in an assemblage

A) Assessing phylogenetic clustering in an assemblage in relation to a species pool

It tests phylogenetic clustering in an assemblage by assessing if its species are more closely related to each other than expected by chance. For this purpose, it calculates the mean phylogenetic distance (MPD) and mean nearest taxon phylogenetic distance (MNTD) and them compares them to MPD/MNTD values for 1,000 randomly generated samples of an equal number of species drawn without replacement from the list of all available species in the species pool. It can be run for multiple trees. The inputs are: (1) phylogenetic tree(s); (2) community matrix; (3) study area mask raster.

Citation and details:

Abellán P, Carrete M, Anadón JD, Cardador L, Tella JL. 2015. Non-random patterns and temporal trends (1912-2012) in the transport, introduction and establishment of exotic birds in Spain and Portugal. Diversity and Distributions, 22: 263-273.

B) Identifying which clades significantly contribute to the phylogenetic structure in an assemblage

For each node in the phylogeny, it is tested whether it has significantly more descendent taxa in a sample than would be expected by chance by means of a randomization test (equivalent to NODESIG function in Phylocom software). Observed patterns are compared to those from random draws of s taxa from the phylogeny terminals where s is the number of taxa in the sample. The inputs are: (1) phylogenetic tree(s); (2) community matrix.

Citation and details:

Abellán P, Carrete M, Anadón JD, Cardador L, Tella JL. 2015. Non-random patterns and temporal trends (1912-2012) in the transport, introduction and establishment of exotic birds in Spain and Portugal. Diversity and Distributions, 22: 263-273.
