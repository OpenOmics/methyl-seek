#!/usr/bin/env RScript
# Title: Tissue of Origin Measurements
# Authors: Alexandre Pellan Cheng
# Brief description: Creates a text file with the tissue of origin measurements (aka mixing parameters (mp)) for
# sample. References are loaded and have already assumed to be prepped (via references.500.rds)

# setwdInitialize -------------------------------------------------------------------------------------------------
rm(list=ls())
# Load libraries ---------------------------------------------------------------------------------------------
library(Matrix)
library(matrixcalc)
library(stats)
library(data.table)
library(limSolve)
library(parallel)

mixing_parameters_function<-function(b, A, other, sample_name, sum_to_one){
  #Accounting for missing tissues
  if(other){
    A$other<-1
    G<-diag(ncol(A))
    g<-rep(-1,ncol(A))
    g[ncol(A)]<-0
    G<-rbind(G, g)
    h<-rep(0,ncol(A))
    h[ncol(A)+1]<-(-1)
    E<-NULL
    f<-NULL
    sol<-lsei(A=A, B=b, G=G, H=h, E=E, F=f, type=2, fulloutput = T) #calls solve.QP from package quadprog
    sol.df<-data.frame(sol$X)
    sol.df$tissue<-rownames(sol.df)
    #sol.df<-sol.df[!rownames(sol.df) %in% "other",]
    if (sum_to_one){
      normalizing_factor<-sum(sol.df$sol.X)
      sol.df$sol.X<-sol.df$sol.X/normalizing_factor
    }
    rel_error<-(sol$solutionNorm/nrow(A))
    rel_error<-data.frame(rel_error, "rel_error")
    colnames(rel_error)<-c(sample_name, "tissue")
  }
  #Not Accounting for missing tissues
  if(!other){
    G<-diag(ncol(A))
    h<-rep(0,ncol(A))
    E<-(rep(1,ncol(A)))
    f<-1
    sol<-lsei(A=A, B=b, G=G, H=h, E=E, F=f, type=1, fulloutput = T)
    sol.df<-as.data.frame(sol$X)
    sol.df$tissue<-rownames(sol.df)
    if (sum_to_one){
      normalizing_factor<-sum(sol.df$sol.X)
      sol.df$sol.X<-sol.df$sol.X/normalizing_factor
    }
    rel_error<-(sol$solutionNorm/nrow(A))
    rel_error<-data.frame(rel_error, "rel_error")
    colnames(rel_error)<-c(sample_name, "tissue")
  }
  colnames(sol.df)[1]<-sample_name
  sol.df<-rbind(sol.df, rel_error)
  return(sol.df)
}

# Source custom function -------------------------------------------------------------------------------------
# Functions for tissue of origin measurement
group_by_celltype <- function(df){
  tmpA <- data.frame(t(df))
  tmpA$tissue <- gsub('_.*', '', colnames(df))
  new_df <- aggregate(.~tissue, tmpA, mean)
  new_names <- new_df$tissue
  new_df$tissue<-NULL
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- new_names
  new_df <- as.matrix(new_df)
  return(new_df)
}

# Main function to call ----------------------------------------------------------------------------
mp_function<-function(reference_and_samples, other, nt, sum_to_one, group_tissues, sample_name){
  A<-as.data.frame(reference_and_samples[,4:(nt+3)])
  if (group_tissues){
    A<-group_by_celltype(A)
  }
  
  b<-reference_and_samples[,(nt+4):ncol(reference_and_samples)]
  percentages <- mixing_parameters_function(b, A, other, sample_name, sum_to_one)
  p <- data.frame(t(percentages))[1,]
  colnames(p)[ncol(p)]<-'rel_error'
  p$sample <- rownames(p)
  return(p)
}


# Source custom function ------------------------------------------------------------------------------------- 
# Functions for tissue of origin measurement
group_by_celltype <- function(df){
  tmpA <- data.frame(t(df))
  tmpA$tissue <- gsub('_.*', '', colnames(df))
  new_df <- aggregate(.~tissue, tmpA, mean)
  new_names <- new_df$tissue
  new_df$tissue<-NULL
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- new_names
  new_df <- as.matrix(new_df)
  return(new_df)
}

# mixing parameters function
mp_function <- function(ref_and_sam, other, sum_to_one, group_tissues, sample_name){
  A <- ref_and_sam[, !(c('chr', 'start', 'end', 'pmeth'))]
  b <- ref_and_sam$pmeth
  
  if (other){
    A$other <- 1
  }
  
  if (group_tissues){
    A <- group_by_celltype(A)
  }else{
    A <- as.matrix(A)
  }
  
  mod <- nnls(A, b)
  percentages <- data.frame(mod$X)
  
  colnames(percentages) <-c('prediction')
  percentages$prediction <- as.numeric(as.character(percentages$prediction))
  
  if (sum_to_one){
    percentages$prediction <- percentages$prediction/sum(percentages$prediction)
  }
  percentages <- data.frame(t(percentages))
  colnames(percentages) <- colnames(A)
  percentages$sample <- sample_name
  return(percentages)
  
}

# Command line arguments -------------------------------------------------------------------------------------

# ARGUMENTS MUST BE SPECIFIED IN THE FOLLOWING ORDER
# 1- sample file path
# 2- reference file path
# 3- output file path
# 4- lookup table (reference methylomes)
# 5- sample name
# 6- sum to one: do we want to sum to one? TRUE or FALSE
# 7- other: do we want the possibility of an error absorbing term? TRUE or FALSE
# 8- group by cell type: do we average the references by group first? TRUE or FALSE
# removals: what tissues in the reference to you want to remove before deconvoluting? no spaces, front slash separated
#               ex: bcell_1/t_cell3

args <- commandArgs(trailingOnly = TRUE)

sam <- args[[1]]
refs <- args[[2]]
out <- args[[3]]
lookup <- args[[4]]
sample_name <- args[[5]]
sum_to_one <- args[[6]]
other <- args[[7]]
GBC <- args[[8]]


if (!(sum_to_one == "TRUE" | sum_to_one == "FALSE")){
  stop('Your sum_to_one parameter is wrong. TRUE or FALSE with caps.')
}else{
  sum_to_one <- as.logical(sum_to_one)
}
if (!(other == "TRUE" | other == "FALSE")){
  stop('Your other parameter is wrong. TRUE or FALSE with caps.')
}else{
  other <- as.logical(other)
}
if (!(GBC == "TRUE" | GBC == "FALSE")){
  stop('Your group_by_celltype parameter is wrong. TRUE or FALSE with caps.')
}else{
  GBC <- as.logical(GBC)
}

if ((length(args)) ==9){
  removals <- args[[9]]
  removals <- gsub('/', '|', removals)
}else{
  removals <- 'thiswillneverbeatissuename' 
}

# Read files and format column names -----------------------------------------------------------------------
sam.df <- fread(sam)
colnames(sam.df)<- c('chr', 'start', 'end', 'pmeth')
sam.df$pmeth <- sam.df$pmeth/100

refs.df <- fread(refs)
coln<-fread(lookup, header=T)
colnames(refs.df)<-c("chr", "start", "end", coln$tissue_name)

refs.df <- data.table(data.frame(refs.df)[, !grepl(removals, colnames(refs.df))])

# Initialize variables ----------------------------------------------------------------------------------------
num.tissues<-ncol(refs.df)-3 # 3 columns are for chromosome, start and end

# Table preparation ------------------------------------------------------------------------------------------
ref_and_sam <- merge(refs.df, sam.df, by=c('chr', 'start', 'end'))
ref_and_sam <- ref_and_sam[complete.cases(ref_and_sam)]

# Measuring tissue of orign ----------------------------------------------------------------------------------
tissues_of_origin <- mp_function(ref_and_sam = ref_and_sam,
                                 other=other,
                                 sum_to_one = sum_to_one,
                                 group_tissues = GBC,
                                 sample_name=sample_name)


#tissues_of_origin
fwrite(tissues_of_origin, file = out, quote=FALSE, sep="\t")

