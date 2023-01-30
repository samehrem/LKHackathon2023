####################Running GWAS on our phenotypes
##Libraries
## load needed packages
library(Matrix)
library(MASS)
library(ggplot2)
library(qqman)
library(ggplot2)
library(cowplot)
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
#BGI_Sat_kinship.out
#phenotypes_for_GWAS.out
#BGI_Sat_kinship.out
#sat.snps.out

#############################

## During GWAS we use linear models to test the association of a genetic variant (SNP,PAV,CNV,kmer) 
## with the trait of interest. For this workshop we run GWAS with SNPs.


##One important step is calculating a kinship matrix. The kinship matrix is used during the GWAS to correct
##for false positive associations by taking into account population structure. We can get an estimate of kinship
##calculating the covariance of the SNPs.


load("./Data/BGI_Sat_kinship.out")
heatmap(letkin)

## We prepared our phenotype input and saved it as phenotypes_for_GWAS.out

## Below you will find the GWAS script we will use. First you run the function as a whole (Ctrl+Enter on start of line 52)

GWAS <- function(genotypes, trait, phenotype.name, kinship, out.dir,
                 give.pval.output.in.R = F, maf.thr = 0.05, snp.info) {
  
  ## Prep trait/phenotypes
  phenotype <- as.character(phenotype.name)
  letkin <- kinship
  usemat <- genotypes
  selc <- !is.na(trait) #Selects the lines with an observation (removes lines that have an NA)
  trait.names <- trait[selc]
  use.trait <- trait[selc]
  print("Traits are selected.")
  print(paste("The phenotype ID is ", phenotype,".", sep=""))

  ## Filter in usemat and kinship object for mapping
  usemat <- usemat[,selc]
  letkin <- letkin[selc,selc]
  
  
  ## Filter again according to MAF we chose
  threshold <- round(ncol(usemat)*(1-maf.thr), digits=0)
  maf.filter.quick <- apply(usemat == 1,1,sum)> threshold  | apply(usemat == 3,1,sum) > threshold
  print(paste0("SNPs falling within MAF > ",(maf.thr)*100,"%" ))
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
  ##Calculate narrow sense heritability
  herit.mod <- lme4qtl::VarProp(mod)
  herit.mod <- herit.mod$prop[1]
  
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
  lod[snp.selc] <- gassoc_gls$tab$pval##Here we add the SNPs we tested, the rest is 0
  
  zscore <- rep(0,nrow(snp.info))
  zscore[snp.selc] <- gassoc_gls$tab$zscore ##Here we add the SNPs we tested, the rest is 0
  mrkno <- which.max(lod)
  
  se <- rep(0,nrow(snp.info))
  se[snp.selc] <- gassoc_gls$tab$se ##Here we add the SNPs we tested, the rest is 0
  
  b <- rep(0,nrow(snp.info))
  b[snp.selc] <- gassoc_gls$tab$b ##Here we add the SNPs we tested, the rest is 0
  
  
  
  ###Save as integers
  gassoc_gls <- snp.info
  gassoc_gls$pval <- lod
  gassoc_gls$zscore <- as.integer(zscore*10000)
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
main.dir <- paste(ext.dir,"GWAS_Results/",sep="")
dir.create(main.dir)

###INPUT

#### Load phenotype data
pheno <- load("./Data/phenotypes_for_GWAS.out")
pheno <- eval(parse(text=pheno))
rm(phenotypes.mean)
base.dir <- paste(main.dir, "BGI_",sep="") ##Indicate which variants were used

#Load genotype object for GWAS mapping
usemat <- load("./Data/sat.snps.out")
usemat <- eval(parse(text=usemat))
snp.info <- usemat[,1:3]
usemat <- data.matrix(usemat[,-c(1:3)])
rm(sat.snps)


#Load kinship matrix
load("./Data/BGI_Sat_kinship.out")

###Input phenotype
trait <- rownames(pheno)[1]
new.dir <- paste(base.dir,trait,sep="")
dir.create(new.dir)
pheno <- pheno[,match(colnames(letkin),colnames(pheno))]
pheno <- pheno[,!is.na(colnames(pheno))]

## Prep sets for GWAS
log.file.w <- file(paste(new.dir,"/","BGI_",trait,"_warning.log",sep=""),open="wt")
sink(file=log.file.w,type="message")

letkin <- letkin[names(pheno[trait,]),names(pheno[trait,])] #In case we do not have information for all lines with this phenotype
usemat<- usemat[,names(pheno[trait,])] #In case we do not have information for all lines with this phenotype


GWAS.output <- GWAS(genotypes = usemat, trait = as.vector(pheno[trait,]), phenotype.name = trait, kinship=letkin, out.dir=new.dir,
                    maf.thr = 0.05,give.pval.output.in.R = T,snp.info=snp.info)

sink(type="message")
close(log.file.w)

print(paste("GWAS finished. Phenotype is ",trait,sep=""))
#print(paste(nrow(pheno) - i," traits of ",nrow(pheno)," to go.",sep=""))
lifecycle::last_lifecycle_warnings()


#Q3: Plot the manhattan plots for the GWAS run 

load(paste("./GWAS_Results/BGI_",trait,"/GWAS_result_",trait,".out",sep=""))
bf <- -log10(0.05/nrow(gassoc_gls)) #Bonferroni threshold
gassoc_gls$POS <- gassoc_gls$POS/1000000  ##We divide here by 1 million to get Mbp instead of bp, that is not needed, I just like the look.
gassoc_gls.topl <- gassoc_gls[-log10(gassoc_gls$pval) >1,]

manhattan(gassoc_gls.topl,snp="SNP",chr="CHR",bp = "POS",p = "pval",logp = T,suggestiveline = F,
          genomewideline = bf,annotatePval = bf,col = c("royalblue4","skyblue"))

###We can also check the narrow sense heritability estimate

load(paste("./GWAS_Results/BGI_",trait,"/Heritability_estimate_",trait,".out",sep=""))
print(paste("h2 is estimated to be: ",round(herit.mod*100,2),"%",sep=""))

##We can also look at how the TopSNP behaves in relation to the phenotype, since we expect a linear relationship
##Between alleles and phenotype.
load(paste("./GWAS_Results/BGI_",trait,"/GWAS_cofac.out",sep=""))

to.pl <- as.data.frame(cbind(pheno[1,match(names(cofac),colnames(pheno))],cofac))
to.pl[,2] <- as.character(to.pl[,2])
p <- ggplot(to.pl, aes(x=cofac, y=V1, fill=cofac)) +
  geom_boxplot(notch=F)+
  geom_jitter(position=position_jitter(0.1))+
  ylab("Phenotype value")+
  xlab("Allele status")+
  ggtitle("Phenotype value against allele status of top SNP")+
  annotate("label",
           x = 1:length(table(to.pl$cofac)),
           y = aggregate(V1 ~ cofac, to.pl, median)[ , 2],
           label = table(to.pl$cofac),
           col = "black",
           fontface="bold")+
  theme(legend.position = "none")

p