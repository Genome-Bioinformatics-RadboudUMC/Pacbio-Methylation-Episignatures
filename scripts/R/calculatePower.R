#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(iterators)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
case_samples = as.numeric(args[1])
control_samples = as.numeric(args[2])

print(case_samples)
print(control_samples)

case_distribution <- fread('cpgSimulationTable_PacBioCases.tsv')
control_distribution <- fread('cpgSimulationTable_PacBioControls.tsv')
epic_array_sites <- fread('Hampstead_KMT2A_MethylationDifferences_withCoordinates_hg38.txt')
epic_array_sites[,id:=paste0(chr,'_',start,'_',end)]

sadikovic_signature_sites <- fread('Sadikovic_KMT2A_signature_hg38.tsv')
hampstead_signature_sites <- fread('Hampstead_KMT2A_consensus_hg38.tsv')
signature_sites <- merge(sadikovic_signature_sites,hampstead_signature_sites,by=c('V1','V2','V3'),all.x=FALSE,all.y=FALSE)
signature_sites[,id:=paste0(V1,'_',V2,'_',V3)]
setkey(signature_sites,id)

case_distribution <- case_distribution[cpg_site %in% epic_array_sites[,id]]
control_distribution <- control_distribution[cpg_site %in% epic_array_sites[,id]]
distribution <- merge(case_distribution,control_distribution,by='cpg_site',all.x=FALSE,all.y=FALSE)
distribution <- distribution[!is.na(sd_meth.x) & !is.na(sd_meth.y)]
setkey(distribution,cpg_site)
distribution[,signature_site:=FALSE][signature_sites,signature_site:=TRUE]

# tests <- data.table(expand.grid(c(10:100),c(10:100)))

readWilcox <- function( case_reads, control_reads, case_mean, case_sd, control_mean, control_sd ) {
  
  case_probabilities <- pmin(pmax(rnorm(case_reads, mean = case_mean, sd = case_sd),0),100)
  control_probabilities <- pmin(pmax(rnorm(control_reads, mean = control_mean, sd = control_sd),0),100)
  return(wilcox.test(case_probabilities,control_probabilities)$p.value)
  
}

# numCores = detectCores()
# registerDoParallel(numCores)
  
distribution[,sampled_cases:=case_samples]
distribution[,sampled_controls:=control_samples]
distribution[,p:=readWilcox(case_samples,control_samples,mean_meth.x,sd_meth.x,mean_meth.y,sd_meth.y),by=cpg_site]
distribution[,corrected_p:=p.adjust(p,method='fdr',n=20000000)]
distribution[,p_pos:=distribution[signature_site==TRUE & corrected_p < 0.05,.N]/distribution[signature_site==TRUE,.N]]
distribution[,p_neg:=distribution[signature_site==FALSE & corrected_p < 0.05,.N]/distribution[signature_site==FALSE,.N]]
distribution[,power:=p_pos - p_neg]
colnames(distribution) <- c('cpg_site','case_mean','case_sd','control_mean','control_sd','signature_site',
                              'sampled_cases','sampled_controls','p','corrected_p','p_pos','p_neg','power')
samplingTable <- unique(distribution[,c('sampled_cases','sampled_controls','p_pos','p_neg','power')])

fwrite(samplingTable,paste0('powerTable_sharedSites_',
       case_samples,'_',control_samples,'.tsv'),
       col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
