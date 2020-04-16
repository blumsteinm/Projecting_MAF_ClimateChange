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
## Note: this is a small subset of the full dataset used to illustrate
##       an example script and thus will not perfectly repicate the 
##       results of the original study
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
##         **Note: we chose to use PCAs of climatic data rather than the raw
##         data due to high colinerity amongst variables. 
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

## v. Make dataframe of results from all PCs for past and future climate for use in the CCA
past_climate_pcs   <- data.frame(Population = past_climate$Population, past_climate_pca$x)
future_climate_pcs <- data.frame(Population = future_climate$Population, future_climate_pca)

## vi. Update the column names for easy identification of climatic variables
colnames(past_climate_pcs)[2:ncol(past_climate_pcs)] <- paste0("Climate_", colnames(past_climate_pcs)[2:ncol(past_climate_pcs)] )
colnames(future_climate_pcs)[2:ncol(future_climate_pcs)] <- paste0("Climate_", colnames(future_climate_pcs)[2:ncol(future_climate_pcs)] )

###############################################################################
## Step 3: Fit a Canonical Correspondence Analysis (CCA) to distributions of 
##         current allele frequencies versus climate normals
###############################################################################

## i. Merge the PC data from the Climate Normals with the allele frequncy data
##    for use in the CCA model
cca_data <- merge(allele_frequencies, past_climate_pcs)

## ii. Find the index of the climate data and the 
climate_index <- grep("Climate", colnames(cca_data))
allele_index  <- grep("X0", colnames(cca_data))

## ii. Run CCA with past climate normals and current allele frequencies by populations.
##     CCAs were initially used in community ecology to predict species distributions
##     against environment, here we swap allele frequncies for speices. 
cca_all  <- cca(cca_data[,allele_index] ~ ., data = cca_data[,climate_index], permutations = 4000)
cca_null <- cca(cca_data[,allele_index] ~ 1, data = cca_data, permutations = 4000)

## iii. Drop environmental predictors that are collinear/non-significant via step-wise running model
revised_cca <- ordistep(cca_null, scope=formula(cca_all), direction = "both")
revised_cca$anova

###############################################################################
## Step 4: Assess model fit and accuracy
##############################################################################

## i. Note that "Total Inertia" is the total variance in allele frequency distributions. "Constrained Inertia" 
##     is the variance explained by the environmental variables (gradients matrix). The "Proportion" values represent the
##     percentages of variance of allele frequency distributions explained by Constrained (environmental) and Unconstrained variables. 
##     Eigenvalues of constrained and unconstrained axes represent the amount of variance explained by each CCA axis 
print("Our CCA model explains...")
cat(round(revised_cca$CCA$tot.chi
          /revised_cca$tot.chi*100), "% of variation in MAF data", "\n")

## ii. Perform significance test to see if our model explains variation in allele frequncies more than expected 
##    by chance (p = 0.05)
anova.cca(revised_cca) ##*gives a pseudo-F via permutation

## iii. Visually compare predicted MAFs to actual MAF. Note in this example
##      some populations are poorly predicted given the limited number of loci
##      used for sample code. If found in real analysis, should consider dropping
##      these populations due to poor predictive ability. 
predict_current <- predict(revised_cca, newdata = future_climate_pcs)

par(mfrow = c(1,1), mar = c(3, 3, 0.1, 0.1), mgp = c(1.5, 0.5, 0))
plot(NULL, xlim = c(0,0.75), ylim = c(0,0.75), ylab = "Predicted Allele Frequency", xlab = "Actual Allele Frequency")
for(k in 1:length(populations)){
  points(predict_current[k,] ~ as.numeric( allele_frequencies[k,2:ncol(allele_frequencies)] ), pch = 16, col = "gray80")
  m <- lm(predict_current[k,] ~ as.numeric( allele_frequencies[k,2:ncol(allele_frequencies)] ))
  abline(m, lwd = 2, col = "gray80")
  text(0.3, 0.75 - (k * 0.02), paste(populations[k], ":", round(summary(m)$r.squared,2)), adj = c(1,0) )
}
abline(a = 0, b = 1, lwd = 2)

###############################################################################
## Step 5: Predict future allele frequencies given climate change using CCA
##         model and calculate summary statistics 
###############################################################################

## i. Use the CCA model of MAF predicted by climate PCs to predict MAF under 
##    future climate 
MAF_future <- predict(revised_cca, newdata = future_climate_pcs)

## ii. Calculate the difference between current and predicted MAFs
MAF_Change <- MAF_future - allele_frequencies[2:ncol(allele_frequencies)]

## iii. Calculate the average shift across all loci by population 
MAF_Change <- rowMeans( abs(MAF_Change) )

## iv. Calculate the number of loci that are missing the minor allele by pop.
MAF_Missing <- rowSums( ifelse(allele_frequencies[2:ncol(allele_frequencies)] < 0.01, 1, 0) ) / (ncol(allele_frequencies) -1) 

## v. Visualize the results 
par(mfrow = c(1,2), mar = c(3, 3, 0.1, 0.1), mgp = c(1.5, 0.5, 0))
plot(MAF_Missing ~ past_climate$Lat, pch = 16, xlab = "Latitude", ylab = "Proportion Loci Missing Minor Allele: Stems", col = "forestgreen")
plot(MAF_Change  ~ past_climate$Lat, pch = 16, xlab = "Latitude", ylab = "Average Prediced MAF Increase: Stems", col = "steelblue")















