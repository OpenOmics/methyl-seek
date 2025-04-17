#!/usr/bin/env Rscript
library(methylKit)
library(dplyr)
library(annotatr)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Rnorvegicus.v79)

args = commandArgs(trailingOnly=TRUE)

datadir <- args[1] ##datadir <- "/data/NHLBI_IDSS/rawdata/NHLBI-6/methylKit/"
methfile <- args[2] ##methfile <- "samples.txt"
outdir <- args[3] ##outdir <- "/data/NHLBI_IDSS/rawdata/NHLBI-6/Amos/"
min_cov <- args[4] ##min_cov <- 5
win_size <- args[5] ##win_size = 1000


##setwd("/data/NHLBI_IDSS/rawdata/NHLBI-6/methylKit")
setwd(datadir)
temp <- read.delim(methfile) ##temp <- read.delim(paste0(datadir, "samples.txt"))
file.list <- as.list(temp$samples)
sample.id = as.list(temp$id)
treatment = c(rep(0, table(temp$group)[1]), rep(1,table(temp$group)[2])) 

myobj=methRead(file.list, sample.id= sample.id, assembly="hg38", treatment=treatment, context="CpG")##This keeps everything regardless of the coverage
filtered.myobj=filterByCoverage(myobj,lo.count=min_cov,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)##removing bases <5x coverage and those with extremely low read coverage (>99.9th percentile) due to PCR bias

##normalize coverage
norm.myobj<- normalizeCoverage(filtered.myobj, method = "median")

# ## -----------------------------------------------------------------------------
# pdf(paste0(outdir, "methQC_",  paste0(min_cov, "X_Coverage"),  ".pdf"), width=9, height=8.5)
# par(mfrow = c(2, 2))
# getMethylationStats(norm.myobj[[1]],plot =TRUE , both.strands = F)
# getCoverageStats(norm.myobj[[1]],plot=TRUE,both.strands=FALSE)
# getMethylationStats(norm.myobj[[length(temp$group)-1]],plot =TRUE , both.strands = F)
# getCoverageStats(norm.myobj[[length(temp$group)-1]],plot=TRUE,both.strands=FALSE)
# dev.off()

##########################################################################
#########Regional analysis########
##########################################################################
tiles = tileMethylCounts(norm.myobj, win.size=win_size,step.size=win_size, cov.bases = 3) 
meth_reg = methylKit::unite(tiles, destrand=TRUE, min.per.group = min(table(temp$group))) 
write.table(meth_reg, file = paste0(outdir, "Region_CT_counts.txt"), sep = "\t", row.names = TRUE, quote = F)

## Sample clustering- PCA
pdf(paste0(outdir, "Region_PCA_Clustering.pdf"), width=9, height=8)
PCASamples(meth_reg)
dev.off()

##Average percent regional methylation#######
pm_reg <- data.frame(percMethylation(meth_reg, rowids = T))
pm_reg$id<- rownames(pm_reg)

library(data.table)
check <- setDT(pm_reg)[, tstrsplit(pm_reg$id, '[.]', type.convert = T)]
colnames(check)<- c("chr", "start", "end")

pm_reg <- cbind(check, pm_reg) %>% dplyr::select(!id) %>% 
  mutate(locnames= paste0("Region", seq(1:nrow(.)))) #%>% 
  ##tibble::column_to_rownames(., var = "locnames")

write.table(pm_reg, file = paste0(outdir, "Region_percent.meth.txt"), sep = "\t", row.names = TRUE, quote = F)

## -----------------------------------------------------------------------------
#####Calculating methdiff at regional level#####
myDiff_reg = calculateDiffMeth(meth_reg, num.cores =8) 
myDiff_reg$locnames <- paste0("Region", seq(1:nrow(myDiff_reg)))

write.table(myDiff_reg, file = paste0(outdir, "Region_methDiff.txt"), sep = "\t", row.names = TRUE, quote = F)


##########################################################################
#########Annotations to CpG features and genic regions########
##########################################################################

complete_annotation<- function(region, model = c("hg38", "hg19", "mm10", "rn6" )){ 
  dm_regions <- region[, 1:ncol(region)] %>% GRanges(.)
  
  if(model == "hg38"){
    annots = c('hg38_cpgs', 'hg38_basicgenes')
    annotations = build_annotations(genome = 'hg38', annotations = annots)
    # Intersect the regions we read in with the annotations
    dm_annotated = annotate_regions(regions = dm_regions, annotations = annotations,
                                    ignore.strand = TRUE, quiet = FALSE) %>% 
      as.data.frame() %>% mutate(annot.id = str_remove(annot.id, "\\:[^:]*$")) %>% 
      mutate(annot.id = ifelse(annot.id=="inter", "interCGI", annot.id)) %>% 
      distinct(locnames, annot.id, .keep_all=TRUE)
    
    ##Annotating to Ensembl ID and distance to TSS
    annoData <- ChIPpeakAnno::toGRanges(EnsDb.Hsapiens.v86, feature="gene") 
    ensembl = readRDS("/data/NHLBI_IDSS/rawdata/NHLBI-6/Amos/ensembl_hg38.rds")
    
  }else if (model == "hg19") {
    annots = c('hg19_cpgs', 'hg19_basicgenes')
    annotations = build_annotations(genome = 'hg19', annotations = annots)
    # Intersect the regions we read in with the annotations
    dm_annotated = annotate_regions(regions = dm_regions, annotations = annotations,
                                    ignore.strand = TRUE, quiet = FALSE) %>% 
      as.data.frame() %>% mutate(annot.id = str_remove(annot.id, "\\:[^:]*$")) %>% 
      mutate(annot.id = ifelse(annot.id=="inter", "interCGI", annot.id)) %>% 
      distinct(locnames, annot.id, .keep_all=TRUE)
    
    ##Annotating to Ensembl ID and distance to TSS
    annoData <- ChIPpeakAnno::toGRanges(EnsDb.Hsapiens.v75, feature="gene")
    ensembl = readRDS("/data/NHLBI_IDSS/rawdata/NHLBI-6/Amos/ensembl_hg19.rds")
    
  }else if(model == "mm10"){
    
    annots = c('mm10_cpgs', 'mm10_basicgenes')
    annotations = build_annotations(genome = 'mm10', annotations = annots)
    # Intersect the regions we read in with the annotations
    dm_annotated = annotate_regions(regions = dm_regions, annotations = annotations,
                                    ignore.strand = TRUE, quiet = FALSE) %>% 
      as.data.frame() %>% mutate(annot.id = str_remove(annot.id, "\\:[^:]*$")) %>% 
      mutate(annot.id = ifelse(annot.id=="inter", "interCGI", annot.id)) %>% 
      distinct(locnames, annot.id, .keep_all=TRUE)
    
    ##Annotating to Ensembl ID and distance to TSS
    annoData <- ChIPpeakAnno::toGRanges(EnsDb.Mmusculus.v79, feature="gene")
    ensembl = readRDS("/data/NHLBI_IDSS/rawdata/NHLBI-6/Amos/ensembl_mm10.rds")
    
  }else{
    annots = c('rn6_cpgs', 'rn6_basicgenes')
    annotations = build_annotations(genome = 'rn6', annotations = annots)
    # Intersect the regions we read in with the annotations
    dm_annotated = annotate_regions(regions = dm_regions, annotations = annotations,
                                    ignore.strand = TRUE, quiet = FALSE) %>% 
      as.data.frame() %>% mutate(annot.id = str_remove(annot.id, "\\:[^:]*$")) %>% 
      mutate(annot.id = ifelse(annot.id=="inter", "interCGI", annot.id)) %>% 
      distinct(locnames, annot.id, .keep_all=TRUE)
    
    ##Annotating to Ensembl ID and distance to TSS
    annoData <- ChIPpeakAnno::toGRanges(EnsDb.Rnorvegicus.v79, feature="gene") 
    ensembl = readRDS("/data/NHLBI_IDSS/rawdata/NHLBI-6/Amos/ensembl_rn6.rds")
  }
  ##Annotating to Genic and CpG features
  anno <- dm_annotated %>% GRanges(.) %>%  
    ChIPpeakAnno::annotatePeakInBatch(., AnnotationData=annoData, maxgap=-1L,
                        output=c("nearestLocation"), multiple=c(TRUE),
                        PeakLocForDistance=c("end"),FeatureLocForDistance=c("TSS"),
                        select=c("all"),ignore.strand=TRUE) %>% as.data.frame(.)
  
  bm <- biomaRt::getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                       filters='ensembl_gene_id', values=anno$feature, mart=ensembl) %>%
    dplyr::rename(feature = ensembl_gene_id) %>% right_join(anno) %>% 
    distinct(locnames, annot.id, feature, .keep_all=TRUE) %>% 
    mutate_if(is.character, list(~na_if(.,""))) %>% 
    mutate(GeneNames = coalesce(external_gene_name, annot.symbol))
  
  return(bm)
}

##Complete annotation (genic, TSS and CpG features)
mydiff <- data.frame(myDiff_reg) %>% left_join(pm_reg[, -c(1:3)]) %>% as.data.frame()
mydiff_annotated <- complete_annotation(region = mydiff, model = "hg38")

write.table(mydiff_annotated, file = paste0(outdir, "Region_methDiff_annotated.txt"), sep = "\t", row.names = TRUE, quote = F)


