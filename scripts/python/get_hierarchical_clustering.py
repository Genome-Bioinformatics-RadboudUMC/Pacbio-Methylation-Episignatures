#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import datetime

###################################
##
## MAIN
##
###################################

## Parameters parser
parser = argparse.ArgumentParser(description="Get dendrogram of hierarchical clustering between multiple PacBio methylation Bed files")
parser.add_argument('-i','--inmatrix', required=True, help='Path to cases controls datamatrix with methylation probabilities (format: TSV)')
parser.add_argument('-ts', '--testings', required=True, help='Path to testings datamatrix with methylation probabilities (format: TSV)')
parser.add_argument('-cs', '--cases', nargs='+', help='Cases samples IDs')
parser.add_argument('-ct', '--controls', nargs='+', help='Controls samples IDs')
parser.add_argument('-g','--graph', default='dendrogram', choices=['clustermap', 'dendrogram'], help='Type of graph to produce, either dendrogram or clustermap')
parser.add_argument('-t','--targets', help='Path to file containing target cpg positions to restrain on (format: TSV)')
parser.add_argument('-o', '--outlabel', required=True, help='Path to output dendrogram')
args = parser.parse_args()

## initialise timer
start = datetime.datetime.now()

cg_dict = {}

## store target cpg positions
reg_df = pd.read_csv(args.targets, sep = "\t", header=None, usecols=[0], names=["cpg"])
print(reg_df)
print("signature dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## load DataMatrix of cases controls methylation probabilities
inmatrix_df = pd.read_csv(args.inmatrix, sep = "\t")
print("cases controls meth proba df loaded")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## load DataMatrix of testings methylation probabilities
testings_df = pd.read_csv(args.testings, sep = "\t")
print("testings meth proba df loaded")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## restrict full datamatrix on positions from signature
inmatrix_df = inmatrix_df.loc[inmatrix_df['cpg'].isin(reg_df['cpg'].tolist())].reset_index(drop=True)
inmatrix_id = [c for c in inmatrix_df.columns.tolist() if c != "cpg"]
inmatrix_df = inmatrix_df[inmatrix_id]
print("full dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## restrict testings datamatrix on positions from signature
testings_df = testings_df.loc[testings_df['cpg'].isin(reg_df['cpg'].tolist())].reset_index(drop=True)
testings_id = [c for c in testings_df.columns.tolist() if c != "cpg"]
testings_df = testings_df[testings_id]
print(testings_df)
print("testings_df dataframe restricted")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

##set a common dataframe with methylation probabilities only
curr_df = pd.concat([inmatrix_df,testings_df],axis=1)
print(curr_df)
print("final dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## rename samples by status before clustering
curr_df = curr_df.rename(columns={spl: "Female_train_" + str(args.cases.index(spl) + 1) for spl in args.cases})
curr_df = curr_df.rename(columns={spl: "Male_train_" + str(args.controls.index(spl) + 1) for spl in args.controls})
curr_df = curr_df.rename(columns={spl: "Female_test_" + str(testings_id.index(spl) + 1) for spl in [testings_id[0], testings_id[1]]})
curr_df = curr_df.rename(columns={spl: "Male_test_1" for spl in [testings_id[2]]})
print(curr_df)
print("changed sample name by status")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

if args.graph == "dendrogram":
    ## produce data for hierarchical clustering
    df_transposed = curr_df.T.values
    df_dist = pdist(df_transposed)
    df_linked = linkage(df_dist, "ward")
    print("Data formated")
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
    ## plot dendrogram
    plt.figure(figsize=(15,14))
    target_cpg = "chr1-22" if args.targets is None else "selected positions"
    plt.title("Hierarchical clustering based on methylation probabilities on " + target_cpg, size=12)
    dn1 = dendrogram(df_linked, labels = curr_df.columns)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.tight_layout()
    title_plot = args.outlabel + '.dendrogram_plot.png'
    plt.savefig(title_plot)
    print("plot done " + title_plot)
    plt.close()
    print("Hierarchical clustering dendrogram done")
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
elif args.graph == "clustermap":
    plt.figure(dpi=300)
    target_cpg = "signature CpGs"
    ##plt.title("Clustermap of methylation probabilities on " + target_cpg, size=14)
    ##sns.clustermap(curr_df, method='ward', yticklabels=False, cmap='coolwarm', xticklabels=1, figsize=(9, 8), cbar_pos=None)
    sns.clustermap(curr_df, method='ward', yticklabels=False, cmap='coolwarm', xticklabels=1, figsize=(9, 8))
    plt.xticks(size=11)
    plt.tight_layout()
    title_plot = args.outlabel + '.colscale.clustermap.png'
    plt.savefig(title_plot)
    print("plot done " + title_plot)
    plt.close()
    print("Clustermap done")
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
