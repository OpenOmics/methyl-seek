#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#library(devtools)
#install_github("juliedwhite/miamiplot", build_vignettes = TRUE)
library(miamiplot)
library(readr)
library(rtracklayer)
library(tidyverse)
library(ggplot2)



##args[1] = path to bsseq folder
##args[2] = p-value/FDR cutoff = 0.05
##args[3] = number of CpGs to highlight
##args[4] = figure name for miami_plot
##args[5] = figure name for violin_plot
##args[6] = figure name for pie_cpg
##args[7] = figure name for pie_genic

comb_file <- list.files(path = args[1], pattern = "./*wtFDR_unfiltered.bed") %>% lapply(., read.table) %>% do.call(rbind, .)
comb_file <- comb_file[-1, ] %>% filter(!str_detect(V1, "chrom")) %>% setNames(c("chr", "start", "end",  "pval", "t")) %>% 
  mutate_at(c(2:5), as.numeric) %>% mutate(FDR = p.adjust(pval, method = "fdr"), rsid = paste0("CpG", seq(1:nrow(.))), 
                                           pos = (start+end)/2, t = 100*t, chr = str_remove(chr, "chr")) %>% 
  select(c(rsid, chr, pos, t, pval, FDR, start, end)) %>%
  mutate(chr = as.numeric(str_replace_all(chr, "X", "23"))) %>% ##chr column must be numeric
  data.frame(.)
  
###Miami plot
miami_plot = ggmiami(data = comb_file,
        split_by = "t", split_at = 0, p = "pval", 
        upper_ylab = "Positive t values",
        lower_ylab = "Negative t values",
        suggestive_line = NULL,
        genome_line = max(comb_file$pval[comb_file$FDR<args[2]]), genome_line_color = "black",
        upper_highlight = comb_file$rsid[comb_file$FDR<args[2] & abs(comb_file$t)>=args[2]], upper_highlight_col = "rsid",
        upper_highlight_color = "red",
        lower_highlight = comb_file$rsid[comb_file$FDR<args[2] & abs(comb_file$t)>=args[2]], lower_highlight_col = "rsid",
        lower_highlight_color = "deepskyblue") 

###Volcano plot
##reference: https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/
v_plot<- ggplot(comb_file) + geom_point(aes(x = t, y = -log10(FDR)), size = 0.7, 
             colour = ifelse(comb_file$FDR<args[2] & comb_file$t>=args[2],"firebrick3",
                             ifelse(comb_file$FDR<args[2] & comb_file$t<= (-1*args[2]), "deepskyblue", "gray50")))+  
  geom_hline(yintercept=-log10(args[2]), color="black", linetype="dashed", size = 0.5)+
  geom_vline(xintercept= c(-20, 20), color="black", linetype="dashed", size = 0.5)+
  xlab("Methylation difference (%)")+ ylab("-log10(FDR)") + theme_bw(base_size = 10) 


###Pie charts for significant hits (FDR<args[2] and abs(meth.diff)>10%)
##CpG features
######setwd("~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/methyl-seek")
cpg.f <- lapply(c("island_annot.bed", "shores_annot.bed", "shelves_annot.bed"), read.table) %>% 
  setNames(c("Island","Shores","Shelves")) %>% sapply(., nrow)
pie_cpg<- pie(cpg.f, labels = names(cpg.f), main = "CpG features")

##Genic regions
genic.f <- lapply(c("5utr_annot.bed", "promoters_annot.bed", "exons_annot.bed","introns_annot.bed", "3utr_annot.bed"), read.table) %>% 
  setNames(c("5UTR", "Promoters","Exons", "Introns", "3UTR")) %>% sapply(., nrow)
pie_genic <- pie(genic.f, labels=names(genic.f), main = "Genic regions")

###save plots
ggsave(plot = miami_plot, filename = args[4], width = 40,height = 10, units = "cm",device = "png")
ggsave(plot = v_plot, filename = args[5], width = 20,height = 20, units = "cm",device = "png")

par(mfrow=c(1,2))
png(plot = pie_cpg, filename = args[6], width = 20,height = 20, units = "cm",device = "png")
png(plot = pie_genic, filename = args[7], width = 20,height = 20, units = "cm",device = "png")
