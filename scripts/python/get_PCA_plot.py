#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import datetime

###################################
##
## MAIN
##
###################################

## Parameters parser
parser = argparse.ArgumentParser(description="Get Principal Component Analysis from cases and controls groups on signature sites")
parser.add_argument('-i','--inmatrix', required=True, help='Path to cases controls datamatrix with methylation probabilities (format: TSV)')
parser.add_argument('-ts', '--testings', required=True, help='Path to testings datamatrix with methylation probabilities (format: TSV)')
parser.add_argument('-cs', '--cases', nargs='+', help='Cases samples IDs')
parser.add_argument('-ct', '--controls', nargs='+', help='Controls samples IDs')
parser.add_argument('-t','--targets', help='Path to file containing target cpg positions to restrain on (format: TSV)')
parser.add_argument('-o', '--outlabel', required=True, help='Path to output plot')
args = parser.parse_args()

## initialise timer
start = datetime.datetime.now()

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
print(inmatrix_df)
print("train dataframe done")
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
print("merged dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## transpose columns and rows
curr_df = curr_df.set_index("cpg").rename_axis(None).T.reset_index(names=["spl_id"])
print(curr_df)
print("transposed dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## set column status
curr_df["status"] = ["status"]*len(curr_df)
## conditional formating of status
train_spl = args.cases + args.controls
test_spl = [s for s in testings_df.columns.tolist() if s not in train_spl]
curr_df.loc[curr_df["spl_id"].isin(args.cases), "status"] = "Females_train"
curr_df.loc[curr_df["spl_id"].isin(args.controls), "status"] = "Males_train"
curr_df.loc[curr_df["spl_id"].isin([test_spl[0]]), "status"] = "Female1_test"
curr_df.loc[curr_df["spl_id"].isin([test_spl[1]]), "status"] = "Female2_test"
curr_df.loc[curr_df["spl_id"].isin([test_spl[2]]), "status"] = "Male_test"

print(curr_df)
print("final dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## perform PCA analysis
features = curr_df[[c for c in curr_df.columns.tolist() if c not in ["spl_id", "status"]]]

## get unique labels
labels = ["Female1_test", "Female2_test", "Male_test", "Females_train", "Males_train"]

## set the model
pca = PCA()

## Fit the model with features and apply the dimensionality reduction on features
feat_r = pca.fit(features).transform(features)

## plot explained_variance_ratio for each component of PCA
plt.figure(figsize=(6,5), dpi=300)
##plt.title("Percentage of variance explained\n" + "by each component of PCA\n", size=12)
plt.plot([n for n in range(1, len(pca.explained_variance_ratio_) + 1)], pca.explained_variance_ratio_, '-o')
plt.xlabel('n_components', size=12)
plt.ylabel('explained_variance_ratio', size=12)
plt.yticks(size=12)
plt.xticks(size=12)
fig_title = args.outlabel + '.PCA_explained_var_ratio.png'
plt.tight_layout()
plt.savefig(fig_title)
plt.close()
print("PCA explained variance ratio plot done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## plot output with labelles categories
plt.figure(figsize=(6,5), dpi=300)
##plt.title("Principal Component Analysis on PacBio LRS cohort\n", size=12)

for c, s_mark, label_n in zip(["green", "#d7369e", "gold", "red", "blue"], ["^", "^", "^", "s", "+"], labels):
    plt.scatter(feat_r[curr_df["status"] == label_n, 0], feat_r[curr_df["status"] == label_n, 1], c=c, marker=s_mark, label=label_n)

##for s_idx, s_i in zip(list(range(len(curr_df))), curr_df["spl_id"].tolist()):
##    plt.annotate("".join([" ", str(s_idx + 1)]), xy=(feat_r[curr_df["spl_id"] == s_i, 0], feat_r[curr_df["spl_id"] == s_i, 1]))

##for s_idx, s_i in zip(list(range(len(args.cases))), curr_df["spl_id"].loc[curr_df["spl_id"].isin(args.cases)].tolist()):
##    plt.annotate("".join([" ", str(s_idx + 1)]), xy=(feat_r[curr_df["spl_id"] == s_i, 0], feat_r[curr_df["spl_id"] == s_i, 1]))

plt.legend()
plt.xlabel('1st eigenvector', size=12)
plt.ylabel('2nd eigenvector', size=12)
plt.yticks(size=12)
plt.xticks(size=12)
fig_title = args.outlabel + '.PCA_2D_plot.png'
plt.tight_layout()
plt.savefig(fig_title)
plt.close()
print("PCA 2D plot done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)
