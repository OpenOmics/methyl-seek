
#### wrapper to prepare sample.txt for DMR analyses
args <- commandArgs(trailingOnly = TRUE)

wd <- as.character(args[1])
group1=as.character(args[2])
group2=as.character(args[3])
sample.file=as.character(args[4])

setwd(wd)

meta = read.delim(sample.file)
colnames(meta) = c("samples", "group")

contrast = as.character(paste(group1, "-vs-", group2, sep=""))

data = meta[which(meta$group %in% c(group1, group2)),]

CpGpath = NULL
for(i in 1:nrow(data)){
  path = paste(wd, "/CpG/", data$samples[i], "/", data$samples[i], ".bismark_bt2_pe.deduplicated.CpG_report.txt.gz", sep="")
  data$comparisons[i] = contrast
  data$path[i] = path
}

write.table(data, paste(wd, "/contrasts.txt", sep=""), sep = "\t", append = FALSE, row.names = FALSE, quote = FALSE)


