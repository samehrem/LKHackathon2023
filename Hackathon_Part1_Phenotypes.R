###Draft Hackathon

#### Part 1: Working with phenotypes
###########LIBRARIES
library(openxlsx)
library(ggplot2)
##########

###LIST OF FILES NEEDED FOR THIS PART
#Hackathon_Phenotypes_Rep1.xlsx
#Hackathon_Phenotypes_Rep2.xlsx
#Species_info.xlsx
############


####Loading in your raw phenotype data from excel files
raw.phenotypes.rep1 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/Hackathon_Phenotypes_Rep1.xlsx")
raw.phenotypes.rep2 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/Hackathon_Phenotypes_Rep2.xlsx")

##Combining replicates
raw.phenotypes <- rbind(raw.phenotypes.rep1,raw.phenotypes.rep2) #rbind is a function that combines two tables by row

#raw.phenotypes.X <- read.xlsx("PathToExcelSheet") ###Example if you want to load you own phenotypes


####Exploring the phenotypes

#Q1: What phenotypes did we load? How many samples do we have?
colnames(raw.phenotypes.rep1) #Different phenotype names
nrow(raw.phenotypes.rep1) #Number of LK accessions

#Q2: What are the distributions of our phenotypes (by horticultural type)?
lk.ID <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/Species_info.xlsx") ###We load in information about the LK accessions
lk.subtype <- lk.ID$Subgroup_SB[match(raw.phenotypes.rep1$accession,lk.ID$LKID)] ###We extract the subtype information by LK ID
table(lk.subtype) #We look at what subtypes are represented in our set
##CODE PLOT ## Here Basten if you have code for that I am happy if you want to insert it here

#####Calculating broad sense heritability (BSH)
## To get an estimate of how much environment and genetics contribute to trait variation we can estimate the
## broad sense heritability of our trait of interest. If we know the approximate BSH, we can interpret GWAS results better.
## To calculate BSH, we need the raw phenotype values, including all replicates.

pheno.res <- lm(raw.phenotypes$green.trimmed_mean_10.250621~raw.phenotypes$accession)
pheno.res.anova <- anova(pheno.res)
BSH <- pheno.res.anova$`Mean Sq`[1]/(pheno.res.anova$`Mean Sq`[1]+pheno.res.anova$`Mean Sq`[2]) ##Make this a function and then apply?
BSH


## Our phenotypes are not ready for GWAS yet. We need to normalise the trait (if needed) and prepare it
## to be used as an input for our GWAS script. For now we will create a phenotype matrix where rows are the phenotypes
## and columns are the LK accessions


#Q3: What type of measurement do we have? It can be counts, ratios, binary/ordinal values or 
## quantitative measurements.We have replicates so we also have to first take the average across the replicates
###INSERT CODE and EXPLANATION FOR NORMALISATION

phenotypes.diff.mean <- (raw.phenotypes.rep1[,-1]+raw.phenotypes.rep2[,-1])/2
rownames(phenotypes.diff.mean) <- raw.phenotypes.rep1$accession
phenotypes.diff.mean <- t(phenotypes.diff.mean) #We move the table around so each row is a phenotype

#Q4: Save phenotype file as input file for GWAS script.

save(phenotypes.diff.mean,file="/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/phenotypes_for_GWAS.out")


