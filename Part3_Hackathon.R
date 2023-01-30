################ GWAS Follow up
##Libraries
library(ggplot2)
library(cowplot)
##############################

###LIST OF FILES NEEDED FOR THIS PART
#20221208_Lactuca_sativa.annotation_overview.tsv
#GWAS result you generated
#sat.snps.out


#############################
##Now that we have a locus (or more) that is associated with our phenotype of interest, we can go further
##by investigating what genes we find within these loci.

##Load gene annotation file
gene.anno <- read.delim("./Data/20221208_Lactuca_sativa.annotation_overview.tsv")

#Load SNP matrix again
usemat <- load("./Data/sat.snps.out")
usemat <- eval(parse(text=usemat))
snp.info <- usemat[,1:3]
usemat <- data.matrix(usemat[,-c(1:3)])
rm(sat.snps)


##Load GWAS result
load("./GWAS_Results/BGI_height.mean_diff/GWAS_result_height.mean_diff.out")
bf <- -log10(0.05/nrow(gassoc_gls))
gassoc_gls_sig <- gassoc_gls[-log10(gassoc_gls$pval) > bf,] ##Only choosing significant SNPs
peak <- do.call(data.frame,aggregate(POS ~ CHR, gassoc_gls_sig, function(x){ c(min(x),max(x))})) #We extract the peaks per chromosome

genes.per.peak <- apply(peak,1,function (x) {
  chromosome <- as.character(x[1])
  start.pos <- x[2]/1e6-0.05
  end.pos <- x[3]/1e6+0.05
  locus.info <- gene.anno[gene.anno$chromosome.number == chromosome & 
                            between(gene.anno$start.sequence/1e6,start.pos,end.pos),]}) ### /1e6 means divided by 1 million
genes.per.peak <- do.call(rbind,genes.per.peak)


###To zoom into the region and determine specific blocks, we choose the Top SNP of that region, and calculate its 
###Correlation to all other SNPs in that window. SNPs correlating with each other usually indicates a form of linkage.

locus <- 7 ##Choosing the peak on Chromosome 5
gassoc_gls_cor <- gassoc_gls_sig[gassoc_gls_sig$CHR ==as.character(locus),] ##Creating the object for the plot
top.snp <- gassoc_gls_cor[which.max(gassoc_gls_cor$pval),1:3]
top.snp.mat <- usemat[which(snp.info$CHR == as.character(top.snp$CHR)&snp.info$POS/1e6 == top.snp$POS/1e6),]

##Here we use a function to get all SNPs within the window of interest (with 50kb up and downstream included)
snp.mat.to.cor <- as.numeric(apply(peak[peak$CHR==locus,],1,function (x) { 
  chromosome <- as.character(x[1])                                       
  start.pos <- x[2]/1e6-0.05
  end.pos <- x[3]/1e6+0.05
  snp.to.cor <- which(snp.info$CHR == chromosome &
                        between(snp.info$POS/1e6,start.pos,end.pos))
  return(snp.to.cor)
}
))

usemat.cor <- usemat[snp.mat.to.cor,] ##Selecting the SNPs from the genotype matrix
snp.info.cor <- snp.info[snp.mat.to.cor,] ###And their positional info

top.snp.cor <- cor(top.snp.mat,t(usemat.cor)) #Correlating the Top SNP to all others
snp.cor.to.pl <- cbind(as.numeric(top.snp.cor),snp.info.cor)

snp.cor.to.pl <- merge(snp.cor.to.pl, gassoc_gls, by=c("CHR","POS"), all.x=TRUE)
colnames(snp.cor.to.pl)[3]<- "correlation.topsnp"

mhpl <- ggplot(snp.cor.to.pl,aes(POS,-log10(pval)))+
  geom_point(aes(colour = correlation.topsnp),alpha=1)+
  scale_colour_gradient2(low = "chocolate", high="darkgreen",mid="skyblue1")+
  xlab("Position (Mbp)") + ylab("-log10(p)") +
  theme_cowplot()+
  geom_hline(yintercept=7.52, linetype='dotted', col = 'red',size=1)+
  
  theme(panel.border = element_rect(colour = "black",linetype = "solid"), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text=element_text(size=15))

mhpl
anno.plot <- ggplot(genes.per.peak[genes.per.peak$chromosome.number == as.character(locus),])+
  geom_segment(aes(x=start.sequence/1e6, xend=stop.sequence/1e6, y=chromosome.number, yend=chromosome.number,color=type), 
               size=5,alpha=0.8)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  xlab("Position (in Mb)") +
  ylab("Annotation")+
  theme(text=element_text(size=15),axis.text.x = element_text(angle = 70, vjust = 0.5, hjust=0.5))

anno.plot
full.locus.plot <- plot_grid(mhpl,anno.plot,ncol = 1,align = "v",rel_heights = c(2,0.5))
full.locus.plot
