####PART 1 ####### Working with phenotypes

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
raw.phenotypes.rep1 <- read.xlsx("./Data/Hackathon_Phenotypes_Rep1.xlsx")
raw.phenotypes.rep2 <- read.xlsx("./Data/Hackathon_Phenotypes_Rep2.xlsx")

##Combining replicates
#rbind is a function that combines two tables vertically (== by rows)
#If you do this with your own data, make sure that columns are the same order
raw.phenotypes <- rbind(raw.phenotypes.rep1,raw.phenotypes.rep2)
                                  
#raw.phenotypes.X <- read.xlsx("PathToExcelSheet") ###Example if you want to load you own phenotypes


####Exploring the phenotypes

#What phenotypes did we load? How many samples do we have?
colnames(raw.phenotypes.rep1) #Different phenotype names
nrow(raw.phenotypes.rep1) #Number of LK accessions

#What are the distributions of our phenotypes (by horticultural type)?
lk.ID <- read.xlsx("./Data/Species_info.xlsx") ###We load in information about the LK accessions
lk.subtype <- lk.ID$Subgroup_SB[match(raw.phenotypes.rep1$accession,lk.ID$LKID)] ###We extract the subtype information by LK ID
table(lk.subtype) #We look at what subtypes are represented in our set

## lets check the phenotypic distribution.
phenotype <-raw.phenotypes$green.trimmed_mean_10.250621 
horti_types <- lk.ID$Subgroup_SB[match(raw.phenotypes$accession,lk.ID$LKID)] 
to.plot <- data.frame(horti_types,phenotype)
ggplot(to.plot,aes(phenotype))+
  geom_histogram(aes(fill=horti_types),bins = 100)+
  theme_light()

ggplot(to.plot,aes(phenotype))+
  geom_histogram(aes(fill=horti_types),bins = 100)+
  facet_grid(horti_types~.)+
  theme_light()


## lets check the phenotypic distribution over the horticultural types.
phenotype <-raw.phenotypes$green.trimmed_mean_10.250621 
horti_types <- lk.ID$Subgroup_SB[match(raw.phenotypes$accession,lk.ID$LKID)] 
to.plot <- data.frame(horti_types,phenotype)
ggplot(to.plot,aes(horti_types,phenotype))+
  geom_boxplot(aes(fill=horti_types))+
  theme_light()


## lets check the phenotypic distribution over the genotypes.
phenotype <-raw.phenotypes$green.trimmed_mean_10.250621 
genotypes <- lk.ID$LKID[match(raw.phenotypes$accession,lk.ID$LKID)] 
to.plot <- data.frame(genotypes,horti_types,phenotype)
mp.per.g <- aggregate(to.plot$phenotype,list(to.plot$genotypes),mean,na.rm=T)
to.plot$genotypes <- factor(to.plot$genotypes,levels = mp.per.g$Group.1[order(mp.per.g$x)])
ggplot(to.plot,aes(genotypes,phenotype))+
  geom_boxplot(aes(fill=horti_types))+
  facet_grid(.~horti_types,scale="free_x",space="free_x")+
  theme_light()



#####Calculating broad sense heritability (BSH)
## To get an estimate of how much environment and genetics contribute to trait variation we can estimate the
## broad sense heritability of our trait of interest. If we know the approximate BSH, we can interpret GWAS results better.
## To calculate BSH, we need the raw phenotype values, including all replicates.

pheno.res <- lm(raw.phenotypes$green.trimmed_mean_10.250621~raw.phenotypes$accession)
pheno.res.anova <- anova(pheno.res)

BSH <- (pheno.res.anova$`Mean Sq`[1]-pheno.res.anova$`Mean Sq`[2])/(pheno.res.anova$`Mean Sq`[1]+pheno.res.anova$`Mean Sq`[2])
print(paste("H2 is estimated to be: ",round(BSH*100,2),"%",sep=""))

## Our phenotypes are not ready for GWAS yet. We need to normalize the trait (if wanted),remove outliers and prepare it
## to be used as an input for our GWAS script. For now we will create a phenotype matrix where rows are the phenotypes
## and columns are the LK accessions


##What type of measurement do we have? It can be counts, ratios, binary/ordinal values or 
## quantitative measurements.We have replicates so we also have to first take the average across the replicates
phenotypes.mean <- (raw.phenotypes.rep1[,-1]+raw.phenotypes.rep2[,-1])/2
rownames(phenotypes.mean) <- raw.phenotypes.rep1$accession
phenotypes.mean <- t(phenotypes.mean)


##We can create histograms per phenotype to check the distribution again
apply(phenotypes.mean,1,hist,breaks=50)

##We might see some non-normal distribution, and outliers. We can remove outliers that deviate more than 2 standard deviations from the mean.
sd.per.pheno <- apply(phenotypes.mean,1,sd)
mean.per.pheno <- apply(phenotypes.mean,1,mean)

for (i in 1:nrow(phenotypes.mean)){
  outliers <- phenotypes.mean[i,] > mean.per.pheno[i]+2*sd.per.pheno[i]|phenotypes.mean[i,] < mean.per.pheno[i]-2*sd.per.pheno[i]
  phenotypes.mean[i,outliers] <- NA
}
apply(phenotypes.mean,1,hist,breaks=50) ##Check how it changed


#Save phenotype file as input file for GWAS script.

save(phenotypes.mean,file="./Data/phenotypes_for_GWAS.out")