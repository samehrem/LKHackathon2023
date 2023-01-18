 ###Draft Hackathon

#### Part 1: Working with phenotypes
###########LIBRARIES
library(openxlsx)
library(ggplot2)
library(lm)
##########

###LIST OF FILES NEEDED FOR THIS PART

#Insert file names
############

####Loading in your raw phenotype data
raw.phenotypes.day1.rep1 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/phe_sat_rep1_1106.xlsx")
raw.phenotypes.day2.rep1 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/phe_sat_rep1_2506.xlsx")
raw.phenotypes.day1.rep2 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/phe_sat_rep2_1106.xlsx")
raw.phenotypes.day2.rep2 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/phe_sat_rep2_2506.xlsx")
##Combining reps
raw.phenotypes.day1 <- rbind(raw.phenotypes.day1.rep1,raw.phenotypes.day1.rep2)
raw.phenotypes.day2 <- rbind(raw.phenotypes.day2.rep1,raw.phenotypes.day2.rep2)

#raw.phenotypes.X <- read.xlsx("PathToExcelSheet") ###Example script for if people want to upload their own phenotypes.


####Exploring the phenotypes

#Q1: What phenotypes did we load? How many samples do we have?
colnames(raw.phenotypes.day1.rep1) #Different phenotype names
nrow(raw.phenotypes.day1.rep1) #Number of LK accessions

#Q2: What are the distributions of our phenotypes (by horticultural type)?
lk.ID <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/Species_info.xlsx") ###We load in information about the LK accessions
lk.subtype <- lk.ID$Subgroup_SB[match(raw.phenotypes.day1.rep1$accession,lk.ID$LKID)] ###We extract the subtype information by LK ID
table(lk.subtype) #We look at what subtypes are represented in our set
##CODE PLOT ## Here Basten if you have code for that I am happy if you want to insert it here

#Q3: What are the mean and median values per LK accession/Horticultural type/replicate?
mean.d1 <- apply(raw.phenotypes.day1[,-c(1,412)],2,mean)
mean.d2 <- apply(raw.phenotypes.day2.rep1[,-c(1,480)],2,mean)

med.d1 <- apply(raw.phenotypes.day1[,-c(1,412)],2,median)
med.d2 <- apply(raw.phenotypes.day2[,-c(1,480)],2,median)


#####Calculating broad sense heritability (BSH)
## To get an estimate of how much environment and genetics contribute to trait variation we can calculate the
## broad sense heritability of our trait of interest. If we know the BSH, we can interpret GWAS results better.
## To calculate BSH, we need the raw phenotype values, including all replicates.

###@Basten, I hope this is the correct way....

day1.res <- lm(raw.phenotypes.day1$height.mean~raw.phenotypes.day1$accession)
day1.res.anova <- anova(day1.res)
BSH <- day1.res.anova$`Mean Sq`[1]/(day1.res.anova$`Mean Sq`[1]+day1.res.anova$`Mean Sq`[2]) ##Make this a function and then apply?
BSH


## Our phenotypes are not ready for GWAS yet. We need to normalise the trait (if needed) and prepare it
## to be used as an input for our GWAS script.

#Q5: What type of measurement do we have? It can be counts, ratios, binary/ordinal values or quantitative measurements.
###INSERT CODE and EXPLANATION FOR NORMALISATION

#Q6: Save phenotype file as input file for GWAS script.

###Insert CODE for saving

############NEW SCRIPT ####################Running GWAS on our phenotypes
##Libraries
## load needed packages
library(Matrix)
library(MASS)
library(ggplot2)
# for Step 1: linear mixed model (no SNPs)
library(lme4qtl) 
#install.packages("devtools")
#devtools::install_github("variani/lme4qtl")
# for Step 2: association tests
#install.packages("devtools")

library(matlm)
#install.packages("devtools")
#devtools::install_github("lmehrem/matlm")
library(wlm)
#devtools::install_github("variani/wlm")

library(tictoc)
library(dplyr)

##############################

###LIST OF FILES NEEDED FOR THIS PART

#############################

## During GWAS we use linear models to test the association of a genetic variant (SNP,PAV,CNV,kmer) 
## with the trait of interest. For this workshop we run GWAS with SNPs.

#Q1: How are SNPs, CNVs and PAVs distributed across the genome?

##INSERT CODE FOR COMPARISON (histogram across genome) (If we can load all types of variants at once, otherwise we just use SNPs)

##One important step is calculating a kinship matrix. The kinship matrix is used during the GWAS to correct
##for false positive associations by taking into account population structure.

#Q2: Calculating kinship using SNPs. Visualising the population structure, including the horticultural types

##INSERT CODE

## We prepared our phenotype input and saved it as [FILENAME]

#Q2: Import phenotype file 
##load(BLA)

## Below you will find the GWAS script we will use.

###INSERT GWAS code

#Q3: Plot the manhattan plots for the GWAS run 

##INSERT MANHATTAN PLOT CODE




##############NEW SCRIPT############### GWAS Follow up

