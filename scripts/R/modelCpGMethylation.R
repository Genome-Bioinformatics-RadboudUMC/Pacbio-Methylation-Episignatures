#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(iterators)
library(doParallel)

## Model methylation
## Calculate mean and standard deviation of methylation level in control samples
## Use these values as model parameters to simulate CpG sites
control_samples <- list.files('data_dir/',
                              pattern='*.combined.bed$',
                              full.names=TRUE)

numCores = detectCores() - 1
registerDoParallel(numCores)

cpgSimulationTable <- foreach( sample = control_samples, .combine = 'rbind', .inorder = FALSE ) %dopar% {
  
  methylation <- fread(sample)
  methylation[,file:=basename(sample)]
  methylation[,cpg_site:=paste0(V1,'_',V2,'_',V3)]
  methylation <- methylation[V1 %in% paste0('chr',1:22),c('cpg_site','V9','file')]
  colnames(methylation) <- c('cpg_site','methylation_level','file')
  return(methylation)
  
}

cpgSimulationTable <- cpgSimulationTable[, .(mean_meth = mean(methylation_level), sd_meth = sd(methylation_level)), by = cpg_site]

fwrite(cpgSimulationTable,
       'cpgSimulationTable_PacBioControls.tsv',
       col.names=TRUE,
       row.names=FALSE,
       quote=FALSE,
       sep='\t')

## Model methylation
## Calculate mean and standard deviation of methylation level in case samples
## Use these values as model parameters to simulate CpG sites
sample_set <- fread('sampleTablePacBio.txt',header=TRUE)
sample_set <- gsub('.haplotagged','',sample_set[type=='Case',sample])
sample_set <- sample_set[sample_set!='DNA18-08043']

case_samples <- c(list.files('case_dir/',
                              pattern='*.combined.bed$',
                              full.names=TRUE),
                  list.files('case_dir2/',
                             pattern='*.combined.bed$',
                             full.names=TRUE))
case_samples <- case_samples[grepl(paste(sample_set,collapse='|'),case_samples)]

cpgSimulationTable <- foreach( sample = case_samples, .combine = 'rbind', .inorder = FALSE ) %dopar% {
  
  methylation <- fread(sample)
  methylation[,file:=basename(sample)]
  methylation[,cpg_site:=paste0(V1,'_',V2,'_',V3)]
  methylation <- methylation[V1 %in% paste0('chr',1:22),c('cpg_site','V9','file')]
  colnames(methylation) <- c('cpg_site','methylation_level','file')
  return(methylation)
  
}

cpgSimulationTable <- cpgSimulationTable[, .(mean_meth = mean(methylation_level), sd_meth = sd(methylation_level)), by = cpg_site]

fwrite(cpgSimulationTable,
       'cpgSimulationTable_PacBioCases.tsv',
       col.names=TRUE,
       row.names=FALSE,
       quote=FALSE,
       sep='\t')