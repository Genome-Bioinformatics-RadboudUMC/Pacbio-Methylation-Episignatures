#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(iterators)
library(doParallel)
library(LaF)

## Set simulation parameters
# differential_methylation = seq(from=0.1,to=0.9,by=0.01)
# signature_size = seq(from=100,to=3000,by=100)
# num_cases = seq(5,100,by=1)
# num_controls = seq(5,100,by=1)
# simulations <- data.table(expand.grid(differential_methylation,signature_size,num_cases,num_controls))
differential_methylation = 0.5
signature_size = 3000
num_cases = 10
num_controls = 10

## Set fixed parameters
num_cpus = 32
hypomethylated_ratio = 0.8
num_QCpass_cpgs = 20000000
alpha = 0.05

## Randomly select case and control samples
header <- as.character(fread('simulatedSamples.tsv',
                             header=FALSE,nrow=1))
num_samples <- length(header) - 1

set.seed(123)
samples <- sample(header[grepl('sample',header)],num_cases+num_controls,replace=FALSE)
case_samples <- samples[1:num_cases]
control_samples <- samples[!samples %in% case_samples]

chunk_size <- ceiling(num_QCpass_cpgs/(num_cpus*4))
chunks <- c(seq(1,num_QCpass_cpgs,chunk_size),num_QCpass_cpgs+1)
print(paste0('Running Wilcoxon tests over ',length(chunks),' chunks with size ',chunk_size))

## Prepare file connection
col_types <- c('character',rep('numeric',num_samples))
col_names <- c('cpg_site',paste0('sample_',1:num_samples))
laf <- laf_open_csv("simulatedSamples.tsv",
                    sep='\t', column_types=col_types, column_names=col_names, skip=1)

## Open a connection to the data, only reading in chunked rows and columns in samples
## Perform significance testing (Wilcoxon rank-sum test)
numCores = detectCores()
registerDoParallel(numCores)

wilcoxTests <- foreach( i = 1:(length(chunks)-1), .combine = 'rbind', .inorder = TRUE, .maxcombine = 4 ) %dopar% {
  
    chunk_start <- chunks[i]
    chunk_end <- chunks[i+1] - 1

    selected_data <- data.table(laf[chunk_start:chunk_end, c('cpg_site',samples)])

    if( chunk_start==1 ) {
      
      # Set signature sets
      set.seed(123)
      signature_sites <- sample(1:nrow(selected_data),signature_size,replace=FALSE)
      selected_data[,signature_site:=FALSE][signature_sites,signature_site:=TRUE]
      selected_data[signature_site==TRUE,direction_of_effect:=sample(c('Hypomethylated','Hypermethylated'),
                                                                        1,
                                                                        replace=FALSE,
                                                                        prob=c(hypomethylated_ratio,1-hypomethylated_ratio)),by=cpg_site]
      
      # Set differential methylation at signature sites
      selected_data[signature_site==TRUE, (case_samples) := lapply(.SD, function(x, direction_of_effect) {
        v <- ifelse(direction_of_effect=='Hypomethylated',
                    x - (x * differential_methylation),
                    x + (x * differential_methylation))
        ifelse(v < 0, 0,
               ifelse(v > 100, 100, v))
      }, direction_of_effect), .SDcols = case_samples]
      
    } else {
      
      selected_data[,signature_site:=FALSE]
      selected_data[,direction_of_effect:=NA]
      
    }
    
    # Test significance
    selected_data[,p:=wilcox.test( as.numeric(.SD[,mget(case_samples)]), as.numeric(.SD[,mget(control_samples)]) )$p.value,by=cpg_site]
    return(selected_data[, .(cpg_site,signature_site,direction_of_effect,p)])
    
}   

close(laf)
fwrite(wilcoxTests,'test.txt',
       col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

## Leave out sample
## Perform Mann-Whitney U test
## FDR correct over the number of tests
## Calculate Pearson correlation case and control scores on left out sample (observed distribution)
## Permute case and control labels of observed distribution
## Calculate Pearson correlation case and control scores (expected distribution)








