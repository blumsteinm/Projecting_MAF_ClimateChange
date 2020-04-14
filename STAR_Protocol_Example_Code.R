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

## Set working directory to cloned folder (github.com/blumsteinm/Projecting_MAF_ClimateChange)
working_directory <- "/Users/meghs/Dropbox/PhD_Dissertation/Code/Paper_1_Allele_Freq/STAR_Protocol_Example/"
setwd(working_directory)

## Install Missing Packages: Bioconductor (R v 3.5 or higher)
# install.packages("BiocManager")
# library(BiocManager)
# install("snpStats")

## Install Missing Packages: CRAN
# install.packages("data.table")

## Import Packages
library(data.table)
library(snpStats)
library(vegan)

## Load Data
snps_of_interest <- read.csv("Example_Data/Example_SNPs.csv")
genome_filename  <- "Example_Data/POTR_SNPs_Subset"
past_climate     <- read.csv("Example_Data/Climate_Normals_1961_1990.csv")
future_climate   <- read.csv("Example_Data/Climate_Projections_Ensemble2080s.csv")
population_info  <- read.csv("Example_Data/Population_Information.csv")

###############################################################################
## Step 1: Calculate the Minor Allele Frequency in each population for each
##         loci of interest in the genome
###############################################################################


######################################################################
## a.) Pull the allele information from the .bed files for each SNP
#####################################################################

## i. Get the id of each snp of interest - in this case it is the chromosome and position along the chromosome
snp_list <- as.character( snps_of_interest$rs ) 


## ii. From the .bed file of all SNPs in the genome, pull out the information for just the SNPs of interest
allele_subset <- read.plink(paste0(genome_filename, ".bed"), 
                            paste0(genome_filename, ".bim"), 
                            paste0(genome_filename, ".fam"),
                            na.strings = c("0", "-9"), sep = "." , 
                            select.subjects = NULL, select.snps = snp_list)

## iii. Pull out the allele information
alleles <- as.data.frame(allele_subset$genotypes)

## iv. Make identifier column "Genotype" from rownames
alleles            <- data.frame(Genotype = row.names(alleles), alleles)
row.names(alleles) <- NULL

######################################################################
## b.) Calculate the Minor Allele Frequency by Population 
#####################################################################

## i. Merge the sample name (here Genotype is the identifier) with the population information (population name and lat/lon)
alleles <- merge(population_info, alleles)

## ii. Create function for calculating allele frequency (1 is AA, 2 is AB, 3 is BB); in the case of our data, Allele 1/A is always
##     the minor allele in the population 
freq  <- function(alleles = NULL){
          converted_alleles <- sapply(as.numeric( alleles ), function(x) ifelse(x == 1, 1, ifelse(x == 2, 0.5, ifelse(x == 3, 0, NA))))
          allele_frequency  <- sum(converted_alleles)/length(converted_alleles)
          return(allele_frequency)
}

## iii. Calculate the minor allele frequency of allele A by population; remove unused identifier columns  
allele_frequencies <- aggregate(. ~ Population, data = alleles, FUN = freq)

## iv. Remove NA columns of unused identifier information 
allele_frequencies <- allele_frequencies[,colSums(is.na(allele_frequencies))<nrow(allele_frequencies)]


###############################################################################
## Step 2: Define the major axes of climate variation using a principle 
##         components analysis (PCA) 
###############################################################################

## i. Ensure that both climate files have the same population information in the same order
populations  <- as.character( past_climate$Population )
future_index <- which( as.character(future_climate$Population) %in% populations )

## ii. Use "prcomp" from the vegan v.2.5-6 pacakge to run PCA on the past climate normals data
past_climate_pca <- prcomp(past_climate[5:ncol(past_climate)], scale = T, center = T) ## only run on the columns with climate information, not identifier columns

## iii. Use predict to precit the PC values of the future climate data ** the rows must be in the same order as past climate
future_climate_pca <- predict(past_climate_pca, newdata = future_climate[future_index, 5:ncol(future_climate)])

## iv. Visualize PCA
biplot(past_climate_pca, col = c("black", "steelblue"), xlabs = as.character( past_climate$Population ), xlim = c(-0.5, 0.8), ylim = c(-0.5, 0.8) ) 
text(future_climate_pca[,1:2], populations,col = "gray80") ## overlay where populations are projected to be in the PC space in 2080






























