###############################################################################
###############################################################################
##
## 04/20/2020
##
## Written by: Meghan Blumstein (blumsteinm@gmail.com); using R v.3.6.0
##
## Example code for projecting allele frequency change at adatptive loci 
## under climate change. For complete details on the use and execution
## of this methodology, see Blumstein et al. (2020). Current Biology. 
##
## Inputs: 
##      1. 30-year normals climate data for locations of interest (.csv file)
##      2. Projected climate data from Ensemble of GCM models for locations of interest (.csv file)
##      3. Genomic Resequencing Data for the genotypes in the study (.bed file)
##      4. A list of the genomic loci of interest - gathered from a GWAS, QTL study, or the literature (.csv file)
##
## Outputs: 
##     1. Current distribution of allele frequencies at loci of interest (.csv file)
##     2. Projected changes to allele frequencies at loci of interest (.csv file)
##     3. Summary statistics for your sample set (.csv file)
##     4. Visualizations of allele frequncy change summary statistics (.shp file)
## 
###############################################################################
###############################################################################
rm(list = ls()) ## Clears environment 
###############################################################################

## Import Packages
require(data.table)


## Roots
working_directory <- "C:/Users/BlumsteinMeghan/Dropbox/"
setwd(working_directory)


###############################################################################
###############################################################################


## Import Packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("snpStats")
library(snpStats)
require(plyr)














































