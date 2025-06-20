#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
from scipy.stats import mannwhitneyu
from scipy.stats import false_discovery_control
import datetime

###################################
##
## MAIN
##
###################################

## Parameters parser
parser = argparse.ArgumentParser(description="Get methylation signature CpGs from comparison of cases to controls")
parser.add_argument('-i','--inmatrix', required=True, help='Cases Controls DataMatrix of methylation probabilities values (format:TSV)')
parser.add_argument('-m','--metrics', required=True, help='Cases Controls DataMatrix of dispersion metrics values (format:TSV)')
parser.add_argument('-cs', '--cases', nargs='+', required=True, help='Cases samples IDs as used for training')
parser.add_argument('-ct', '--controls', nargs='+', required=True, help='Controls samples IDs used for training')
parser.add_argument('-dmin','--diffmin', type=float, help='Minimal required difference in median methylation percentage between cases and controls groups')
parser.add_argument('-dmax','--diffmax', type=float, help='Maximal required difference in median methylation percentage between cases and controls groups')
parser.add_argument('-p','--pval', required=True, type=float, help='Threshold for adjusted p-values to keep differentially methylated CpGs')
parser.add_argument('-sp','--spread', type=float, help='Set a maximal threshold for standard deviation in cases and controls groups')
parser.add_argument('-om', '--outmatrix', required=True, help='Path to output file with matrix of median meth differences and pvalues')
parser.add_argument('-os', '--outsign', required=True, help='Path to output file with CpG IDs significant in methylation signature of cases')
args = parser.parse_args()

## initialise timer
start = datetime.datetime.now()

## load DataMatrix of methylation probabilities
inmatrix_df = pd.read_csv(args.inmatrix, sep = "\t")
print(inmatrix_df)
print("meth proba df loaded")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## load DataMatrix of dispersion metrics for cases and controls groups
metrics_df = pd.read_csv(args.metrics, sep = "\t")
print(metrics_df)
print("metrics df loaded")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

if args.diffmin is not None:
    ## restrict to CpG sites with median methylation difference between cases and controls groups higher to specified min percentage
    metrics_df = metrics_df.loc[abs(metrics_df["medians_diff"]) >= args.diffmin]
    print(metrics_df)
    print("restrict medians diffmin")
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)

if args.diffmax is not None:
    ## restrict to CpG sites with median methylation difference between cases and controls groups under specified max percentage
    metrics_df = metrics_df.loc[abs(metrics_df["medians_diff"]) <= args.diffmax]
    print(metrics_df)
    print("restrict medians diffmax")
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)

if args.spread is not None:
    ## restrict to CpG sites with standard deviation by group under cutoff
    metrics_df = metrics_df.loc[metrics_df['SD_cases'] <= args.spread]
    print(metrics_df)
    metrics_df = metrics_df.loc[metrics_df['SD_controls'] <= args.spread]
    print(metrics_df)
    print("restrict under SD threshold in both groups")
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## filter the data matrix of methylation probabilities on selected cpgs based on dispersion filters
inmatrix_df = inmatrix_df.loc[inmatrix_df['cpg'].isin(metrics_df['cpg'].tolist())].reset_index(drop=True)
print(inmatrix_df)

## perform Mann-Whitney U rank test between cases and controls methylation probabilities for each CpG site
inmatrix_df["pvalue"] = mannwhitneyu(inmatrix_df[args.cases], inmatrix_df[args.controls], axis=1).pvalue
print(inmatrix_df)
print("finished statistical test")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## adjust p-values to control the false discovery rate
inmatrix_df["adjusted_pvalue"] = false_discovery_control(inmatrix_df["pvalue"])
print(inmatrix_df)
print("p values adjustment done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## keep only CpG sites with adjusted p-values fitting the specified threshold
inmatrix_df = inmatrix_df.loc[inmatrix_df["adjusted_pvalue"] < args.pval]
print(inmatrix_df)
print("filtered on adjusted p values")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## write out short matrix
inmatrix_df["median_cases"] = inmatrix_df[args.cases].median(axis=1)
inmatrix_df["median_controls"] = inmatrix_df[args.controls].median(axis=1)
inmatrix_df["medians_diff"] = inmatrix_df["median_cases"] - inmatrix_df["median_controls"]
inmatrix_df['SD_cases'] = inmatrix_df[args.cases].std(axis=1)
inmatrix_df['SD_controls'] = inmatrix_df[args.controls].std(axis=1)
print(inmatrix_df)
inmatrix_df[["cpg", "median_cases", "median_controls", "medians_diff", "pvalue", "adjusted_pvalue"]].to_csv(args.outmatrix, sep='\t', mode="w", index=False)

## save selected CpG sites in output file
inmatrix_df["cpg"].to_csv(args.outsign, sep='\t', mode="w", header=False, index=False)
