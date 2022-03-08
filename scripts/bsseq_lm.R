#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(bsseq)
library(tidyverse)
library(dplyr)
library(broom)
library(tidyr)
library(parallel)
library(BiocParallel)

bis_files = read.delim(args[2])

bsList = list()
count = 1
for (i in bis_files$path){
  filename=i
  samplename=bis_files$sample[count]
  print(filename)
  print(samplename)
  tt = read.bismark(files=filename , colData = DataFrame(row.names = samplename), strandCollapse=T, rmZeroCov=T, verbose=T)
  ttChr = subsetByOverlaps(tt, GRanges(seqnames=args[1], ranges=IRanges(start=1, end=200*10^7)))
  bsList <- c(bsList, ttChr)
  #assign(samplename, bsList[[length(bsList)]])
  count = count + 1
}

chr1 = combineList(bsList)

########optional to check that PhenoFile samples are ordered same as the bsList object
new<-as.data.frame(colnames(chr1))
names(new)<-"V1"
new2<-gsub(".*/","",new$V1)
new3<-as.data.frame(gsub("\\..*","",new2))
names(new3)<-"V2"
pheno = bis_files[,-4]
pheno=pheno[match(new3$V2, pheno$sample),]
samplesTomatch=pheno$sample
samplesTomatch
new3$V2
table(samplesTomatch==new3$V2)
########optional to check that PhenoFile samples are ordered same as the bsList object ENDS!!!!
pheno = pheno[,c(1,2)]

###assign phenotype 0 and 1

base_pheno = names(table(pheno$group))[1]
pheno$Dx<-ifelse(pheno$group==base_pheno,1,0)
#pheno = pheno % % filter(`1_Covid_vs_Healthy` != "NA")
pheno
pheno = pheno[,c(1,3)]
samples = pheno$sample

pheno2<-as.data.frame(pheno[,2])
head(pheno2)

###exclude low coverage
cov_frac = as.numeric(as.character(args[4]))
cov_level = as.numeric(as.character(args[5]))
sample_min = min(table(pheno$Dx))
cov_cutoff = ceiling(sample_min*cov_frac)

chr1_cov = getCoverage(chr1, type="Cov", what="perBase")
print("chr1_cov done")
keepLoci.ex = which(rowSums(chr1_cov  >= cov_level) >= cov_cutoff)
print("keepLoci done") 
chr1_filtered = chr1[keepLoci.ex,]
print("chr1_filtered done") 

chr1_betavals_df = data.frame(getMeth(chr1_filtered, type="raw"))
print("chr1_betavals done") 
chr1_granges = granges(chr1_filtered)
chr1_granges_df = data.frame(chromosome=seqnames(chr1_granges),startpos=start(chr1_granges),endpos=end(chr1_granges),strand=strand(chr1_granges))
chr1_betavals_granges_df = data.frame(chr1_granges_df, chr1_betavals_df, row.names=paste(chr1_granges_df$chromosome,chr1_granges_df$startpos,sep="_"), check.names=F)

chr1_betavals_granges_df2 = na.omit(chr1_betavals_granges_df)

keep <- apply(chr1_betavals_granges_df2[,c(-1,-2,-3,-4)], 1, function(x) length(unique(x[!is.na(x)])) != 1)
chr1_betavals_granges_df3 = chr1_betavals_granges_df2[keep, ]

print(dim(chr1_betavals_granges_df3))

new<-data.frame(t(chr1_betavals_granges_df3[,5:dim(chr1_betavals_granges_df3)[2]]),pheno[,-c(1)])
ncolu<-dim(new)[2]
names(new)[ncolu]<-"Dx"
chr1_betavals_granges_df2 = as_tibble(new, rownames="Library_ID")
###run lm
chr1_betavals_granges_df3 = pivot_longer(chr1_betavals_granges_df2, cols=starts_with("chr"), names_to="CpGID_chr_bp", values_to="methylation")

chr1_betavals_granges_df4 = chr1_betavals_granges_df3 %>% group_by(CpGID_chr_bp) %>% nest()

print("create dataframe")
head(chr1_betavals_granges_df4)

mod_fun <- function(df) lm(methylation ~ factor(Dx), data=df)

chr1_lm_results = chr1_betavals_granges_df4 %>% mutate(model = map(data,mod_fun))

chr1_lm_results2 = chr1_lm_results %>% group_by(CpGID_chr_bp) %>% transmute(beta = map_dfr(model,tidy)$estimate[2], pvalue = map_dfr(model,tidy)$p.value[2])

print("model fit complete")

chr1_lm_results3 = as.data.frame(do.call(rbind, strsplit(chr1_lm_results2$CpGID_chr_bp,split = "_")))
chr1_lm_results3$V2 = as.numeric(as.character(chr1_lm_results3$V2))
chr1_lm_results3$V3 = chr1_lm_results3$V2 + 1
chr1_lm_results3$V4 = chr1_lm_results2$pvalue
chr1_lm_results3$V5 = chr1_lm_results2$beta
colnames(chr1_lm_results3) = c("chr","start","end","p","beta")

chr1_lm_results4 = chr1_lm_results3[chr1_lm_results3$p != "NaN",]

write.table(chr1_lm_results4, file=args[3], sep = "\t", row.names=FALSE,quote=F)
