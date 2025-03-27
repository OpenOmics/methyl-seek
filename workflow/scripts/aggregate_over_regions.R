#!/usr/bin/env RScript
# Title: aggregate_over_regions.R
# Authors: Tom Hill
# Brief description: merges retained CpG regions into windows

library(data.table)

custom_aggregate <- function(coll){
  s = sum(as.numeric(as.character(coll[!grepl('-', coll)])))
  
  ratio_known <- length(coll[!grepl('-', coll)]) / length(coll)
  
  if (ratio_known > 0){
    frac <- s / length(coll[!grepl('-', coll)])
  }
  else{
    frac <- NaN
  }
  return(frac)
}

args = commandArgs(trailingOnly=TRUE)

input = args[1]
output = args[2]

sbp_MM = fread(input)

sbp_MM = sbp_MM[, c(1,2,3, 7:ncol(sbp_MM)), with=FALSE]

agg_MM = sbp_MM[, lapply(.SD, custom_aggregate), by=c("V1", "V2", "V3")]

#agg_MM = aggregate(.~V1+V2+V3, sbp_MM, custom_aggregate)
agg_MM = agg_MM[complete.cases(agg_MM), ]
fwrite(x = agg_MM, file = output, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
