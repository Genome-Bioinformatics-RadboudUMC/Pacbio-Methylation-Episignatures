#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(iterators)
library(doParallel)

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
num_samples = args[1]
num_QCpass_cpgs = args[2]

# Get CpG sites to simulate
set.seed(123)
print(paste0('Simulating ',num_samples,' samples at ',num_QCpass_cpgs,' QC-pass CpG sites'))
cpgSimulationTable <- fread('simulate_methylation_signatures/cpgSimulationTable.tsv')
cpgSimulationTable <- cpgSimulationTable[!is.na(sd_meth)]
cpgSimulationTable <- cpgSimulationTable[sample(.N,num_QCpass_cpgs,replace=FALSE)]

## Simulate samples
numCores <- detectCores() - 1
registerDoParallel(numCores)
samples <- paste0('sample_',1:num_samples)

simulatedSamples <- foreach( s = samples, .combine = 'cbind', .inorder = FALSE, .maxcombine = 10 ) %dopar% {
  
  cpgSimulationTable[,sample:=rnorm(1,mean=mean_meth,sd=sd_meth),by=cpg_site]
  cpgSimulationTable[,sample:=ifelse(sample > 100, 100,
                                     ifelse(sample < 0, 0, sample))]
  simulatedSample <- data.table(cpgSimulationTable[,sample])
  colnames(simulatedSample) <- s
  return(simulatedSample)
  
}

simulatedSamples <- cbind(cpgSimulationTable[,cpg_site],simulatedSamples)
colnames(simulatedSamples) <- c('cpg_site',samples)

fwrite(simulatedSamples,'simulatedSamples.tsv',
       col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
