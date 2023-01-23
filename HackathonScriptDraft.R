 ###Draft Hackathon

#### Part 1: Working with phenotypes
###########LIBRARIES
library(openxlsx)
library(ggplot2)
##########

###LIST OF FILES NEEDED FOR THIS PART

#Insert file names
############

####Loading in your raw phenotype data
raw.phenotypes.rep1 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/Hackathon_Phenotypes_Rep1.xlsx")
raw.phenotypes.rep2 <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/Hackathon_Phenotypes_Rep2.xlsx")

##Combining replicates
raw.phenotypes <- rbind(raw.phenotypes.rep1,raw.phenotypes.rep2)

#raw.phenotypes.X <- read.xlsx("PathToExcelSheet") ###Example script for if people want to upload their own phenotypes.


####Exploring the phenotypes

#Q1: What phenotypes did we load? How many samples do we have?
colnames(raw.phenotypes.rep1) #Different phenotype names
nrow(raw.phenotypes.rep1) #Number of LK accessions

#Q2: What are the distributions of our phenotypes (by horticultural type)?
lk.ID <- read.xlsx("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/Species_info.xlsx") ###We load in information about the LK accessions
lk.subtype <- lk.ID$Subgroup_SB[match(raw.phenotypes.rep1$accession,lk.ID$LKID)] ###We extract the subtype information by LK ID
table(lk.subtype) #We look at what subtypes are represented in our set
##CODE PLOT ## Here Basten if you have code for that I am happy if you want to insert it here

#Q3: What are the mean and median values per LK accession/Horticultural type/replicate?
##TODO


#####Calculating broad sense heritability (BSH)
## To get an estimate of how much environment and genetics contribute to trait variation we can calculate the
## broad sense heritability of our trait of interest. If we know the BSH, we can interpret GWAS results better.
## To calculate BSH, we need the raw phenotype values, including all replicates.

###@Basten, I hope this is the correct way....

pheno.res <- lm(raw.phenotypes$relblue.mean~raw.phenotypes$accession)
pheno.res.anova <- anova(pheno.res)
BSH <- pheno.res.anova$`Mean Sq`[1]/(pheno.res.anova$`Mean Sq`[1]+pheno.res.anova$`Mean Sq`[2]) ##Make this a function and then apply?
BSH


## Our phenotypes are not ready for GWAS yet. We need to normalise the trait (if needed) and prepare it
## to be used as an input for our GWAS script. For now we will create a phenotype matrix where rows are the phenotypes
## and columns are the LK accessions


#Q5: What type of measurement do we have? It can be counts, ratios, binary/ordinal values or 
## quantitative measurements.We have replicates so we also have to first take the average across the replicates
###INSERT CODE and EXPLANATION FOR NORMALISATION

##Here i have some questions: for height, red/blue values etc we do log2? 
phenotypes.diff.mean <- (raw.phenotypes.rep1[,2:5]+raw.phenotypes.rep2[,2:5])/2
rownames(phenotypes.diff.mean) <- raw.phenotypes.rep1$accession
phenotypes.diff.mean <- t(phenotypes.diff.mean)
#Q6: Save phenotype file as input file for GWAS script.

save(phenotypes.diff.mean,file="/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/phenotypes_for_GWAS.out")

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


##One important step is calculating a kinship matrix. The kinship matrix is used during the GWAS to correct
##for false positive associations by taking into account population structure. We can get an estimate of kinship
##calculating the covariance of the SNPs.

##NOT POSSIBLE ON A LAPTOP so they will load the kinship

load("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/BGI_Sat_kinship.out")
heatmap(letkin)

## We prepared our phenotype input and saved it as [FILENAME]

## Below you will find the GWAS script we will use.

GWAS <- function(genotypes, trait, phenotype.name, kinship, out.dir,
                 give.pval.output.in.R = F, maf.thr = 0.95, snp.info) {
  
  ## Prep trait/phenotypes
  phenotype <- toString(phenotype.name)
  letkin <- kinship
  usemat <- genotypes
  selc <- !is.na(trait) #Selects the lines with an observation (removes lines that have an NA)
  trait.names <- trait[selc]
  use.trait <- trait[selc]
  print("Traits are selected.")
  print(paste("The phenotype ID is ", phenotype,".", sep=""))
  print(max(usemat,na.rm = T))
  
  ## Filter in usemat and kinship object for mapping
  usemat <- usemat[,selc]
  letkin <- letkin[selc,selc]
  print(dim(letkin))
  

  ## Filter again MAF 5%
  threshold <- round(ncol(usemat)*maf.thr, digits=0)
  print(threshold)
  maf.filter.quick <- apply(usemat == 1,1,sum)> threshold  | apply(usemat == 3,1,sum) > threshold ### <-- this is quicker!
  print(paste0("SNPs falling within MAF >= ",(1-maf.thr)*100,"%" ))
  print(table(maf.filter.quick))
  usemat <- usemat[!maf.filter.quick,]
  snp.info <- snp.info[!maf.filter.quick,]
  print("Genotype matrix filtered and transformed.")
  
  ## Prune SNP set
  phe.snp.cor <- cor(use.trait,t(usemat),use = "pairwise") ###Calculate correlation of SNPs
  print("SNP correlation calculated.")
  phe.snp.cor[is.na(phe.snp.cor)] <- 0 ##Set NAs to 0
  
  snp.selc <- abs(phe.snp.cor)>0.3 & !is.na(phe.snp.cor) 
  usemat.pruned <- usemat[snp.selc,] ##Remove SNPs with an absolute correlation lower than 0.3
  print("SNPs pruned")
 

  ### start mapping by making decomposition
  ID <- rownames(letkin) ; length(ID)
  cbind(ID,use.trait)
  mod <- lme4qtl::relmatLmer(use.trait ~ (1|ID), relmat = list(ID = letkin))
  ##Calculate heritability
  herit.mod <- lme4qtl::VarProp(mod)
  V <- lme4qtl::varcov(mod, idvar = "ID")
  V_thr <- V
  V_thr[abs(V) < 1e-10] <- 0
  decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
  W <- decomp$transform
  print("Decomposition of covariance matrix was performed.")
  
  ## make data object for mapping without any extra factors
  nnn <- rep(1,ncol(letkin))
  # class(use.trait)
  # class(nnn)
  # class(usemat)
  ### GWAS with kinship
  gassoc_gls <- matlm(as.numeric(use.trait) ~ nnn, nnn, pred =  t(usemat.pruned), ids = rownames(W), transform = W, batch_size = 4000, verbose = 2,cores = 1,stats_full = T)
  


  ###Add SNPs we didnt test back
  
  lod <- rep(0.9,nrow(snp.info))
  lod[snp.selc] <- gassoc_gls$tab$pval##Here we add the SNPs we tested, teh rest is 0
  
  zscore <- rep(0,nrow(snp.info))
  print(table(is.na(zscore)))
  zscore[snp.selc] <- gassoc_gls$tab$zscore ##Here we add the SNPs we tested, the rest is 0
  print(table(is.na(zscore)))
  mrkno <- which.max(lod)
  
  se <- rep(0,nrow(snp.info))
  se[snp.selc] <- gassoc_gls$tab$se ##Here we add the SNPs we tested, the rest is 0
  
  b <- rep(0,nrow(snp.info))
  b[snp.selc] <- gassoc_gls$tab$b ##Here we add the SNPs we tested, the rest is 0

  
  
  ###Save as integers
  gassoc_gls <- snp.info
  gassoc_gls$pval <- lod
  gassoc_gls$zscore <- as.integer(zscore*10000)
  print(table(is.na(gassoc_gls$zscore)))
  gassoc_gls$se <- as.integer(se*10000)
  gassoc_gls$b <- as.integer(b *10000)
  gassoc_gls$SNP <- paste(gassoc_gls$CHR,gassoc_gls$POS,sep="_")
  save(gassoc_gls,file=paste(out.dir,"/GWAS_result_",phenotype,".out",sep=""))
  save(herit.mod,file=paste(out.dir,"/Heritability_estimate_",phenotype,".out",sep=""))
  cofac <- usemat[mrkno,]
  save(cofac,file=paste(out.dir,"/GWAS_cofac.out",sep=""))
  print("Results saved.")
  if( give.pval.output.in.R ){
    return(gassoc_gls)
  }
}

########################START SCRIPT#########################
ext.dir <- "./"
main.dir <- paste(ext.dir,"/GWAS_Results/",sep="")
dir.create(main.dir)

###INPUT

#### Load phenotype data
pheno <- load("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/phenotypes_for_GWAS.out")
pheno <- eval(parse(text=pheno))
rm(phenotypes.diff.mean)
base.dir <- paste(main.dir, "BGI_",sep="") ##Indicate which variants were used

#Load genotype object for GWAS mapping
usemat <- load("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/sat.snps.out")
usemat <- eval(parse(text=usemat))
snp.info <- usemat[,1:3]
usemat <- data.matrix(usemat[,-c(1:3)])
rm(sat.snps)


#Load kinship matrix
load("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/BGI_Sat_kinship.out")

###Input phenotype
trait <- rownames(pheno)[2] 
new.dir <- paste(base.dir,trait,sep="")
dir.create(new.dir)
print(colnames(pheno))
pheno <- pheno[,colnames(pheno) %in% colnames(letkin)]
print(ncol(pheno))

## Prep sets for GWAS
log.file.w <- file(paste(new.dir,"/","BGI_",trait,"_warning.log",sep=""),open="wt")
sink(file=log.file.w,type="message")

letkin <- letkin[names(pheno[2,]),names(pheno[2,])] #In case we do not have information for all lines with this phenotype
usemat<- usemat[,names(pheno[2,])] #In case we do not have information for all lines with this phenotype

#GWAS(genotypes = usemat_in, trait = as.vector(pheno[trait,]), phenotype.name = trait, kinship=letkin_in, out.dir=new.dir,
     #maf.thr = 0.95,give.pval.output.in.R = F)
# or 
GWAS.output <- GWAS(genotypes = usemat, trait = as.vector(pheno[trait,]), phenotype.name = trait, kinship=letkin, out.dir=new.dir,
                    maf.thr = 0.95,give.pval.output.in.R = T,snp.info=snp.info)

sink(type="message")
close(log.file.w)

print(paste("GWAS finished. Phenotype is ",trait,sep=""))
#print(paste(nrow(pheno) - i," traits of ",nrow(pheno)," to go.",sep=""))
lifecycle::last_lifecycle_warnings()


#Q3: Plot the manhattan plots for the GWAS run 

##
##INSERT MANHATTAN PLOT CODE

library(qqman)
library(dplyr)
load("/Users/6186130/Documents/LettuceKnow/LKHackathon2023/GWAS_Results/BGI_relgreen.mean/GWAS_result_relgreen.mean.out")
bf <- -log10(0.05/nrow(gassoc_gls)) #Bonferroni threshold
gassoc_gls$POS <- gassoc_gls$POS/1000000
gassoc_gls <- gassoc_gls[-log10(gassoc_gls$pval) >1,]

manhattan(gassoc_gls,snp="SNP",chr="CHR",bp = "POS",p = "pval",logp = T,suggestiveline = F,
          genomewideline = bf,annotatePval = bf,col = c("royalblue4","skyblue"))

##############NEW SCRIPT############### GWAS Follow up

##Now that we have a locus (or more) that is associated with our phenotype of interest, we can go further
##by investigating what genes we find within these loci.



##Load gene annotation file
gene.anno <- read.delim("/Users/6186130/Documents/LettuceKnow/Hackathon_2023/Data/20221208_Lactuca_sativa.annotation_overview.tsv")

gassoc_gls_sig <- gassoc_gls[-log10(gassoc_gls$pval) > bf,]
peak <- do.call(data.frame,aggregate(POS ~ CHR, gassoc_gls_sig, function(x){ c(min(x),max(x))}))

test <- apply(peak,1,function (x) {
  chromosome <- as.character(x[1])
  start.pos <- x[2]-0.1
  end.pos <- x[3]+0.1
  locus.info <- gene.anno[gene.anno$chromosome.number == chromosome & 
              between(gene.anno$start.sequence/1e6,start.pos,end.pos),]})
test <- do.call(rbind,test)
