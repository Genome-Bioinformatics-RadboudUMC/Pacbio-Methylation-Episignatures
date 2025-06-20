#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(iterators)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
imprinted_regions_path = args[1]
num_controls = args[2]
num_resamplings = args[3]

control_samples <- list.files('control_data',
                              pattern='*.combined.bed$',
                              full.names=TRUE)

numCores = detectCores() - 1
registerDoParallel(numCores)

imprinted_regions <- fread(imprinted_regions_path,header=FALSE)
setkey(imprinted_regions,V1,V2,V3)

combineFunction <- function(...) {
  return(merge(...,by='cpg_site',all.x=FALSE,all.y=FALSE))
}

set.seed(1234)
simulationTable <- data.table('cpg_site' = character())

simulationTable <- foreach( n = 1:num_resamplings, .combine = combineFunction, .inorder = FALSE, .multicombine = FALSE ) %dopar% {
  
  sampling_controls <- sample(control_samples,num_controls,replace=FALSE)
  
  cpgTable <- data.table('cpg_site' = character())
  
  cpgTable <- foreach( s = sampling_controls, .combine = combineFunction ) %do% {
    
    methylation <- fread(s)
    setkey(methylation,V1,V2,V3)
    overlaps <- foverlaps(methylation,imprinted_regions,nomatch=NULL)
    overlaps[,cpg_site:=paste0(V1,'_',i.V2,'_',i.V3)]
    overlaps <- overlaps[V1 %in% paste0('chr',1:22),c('cpg_site','V9')]
    overlaps[,V9:=V9/100]
    colnames(overlaps) <- c('cpg_site',basename(s))
    return(overlaps)
    
  }

  cpgTable <- cpgTable[, .(sd_meth = sd(.SD)), .SDcols = grep('cpg_site',colnames(cpgTable),invert=TRUE,value=TRUE), by = cpg_site]
  colnames(cpgTable) <- c('cpg_site',paste0('sampling_',n))
  return(cpgTable)
  
}

# simulationTable[,min_sd:=min(.SD), .SDcols = grep('cpg_site',colnames(cpgTable),invert=TRUE,value=TRUE), by = cpg_site]
# simulationTable <- unique(simulationTable[,c('cpg_site','min_sd')])

fwrite(simulationTable,
       paste0('cpgSimulationTable_PacBioControls_downsampled_',num_controls,'.tsv'),
       col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
