#!/usr/bin/python3

import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import datetime


###################################
##
## MAIN
##
###################################

## Parameters parser
parser = argparse.ArgumentParser(description="Build pearson correlation between test samples and cases controls means")
parser.add_argument('-i','--inmatrix', required=True, help='Path to cases controls datamatrix with methylation probabilities (format: TSV)')
parser.add_argument('-ts', '--testings', required=True, help='Path to testings datamatrix with methylation probabilities (format: TSV)')
parser.add_argument('-s','--sign', required=True, help='Path to the file containing signature sites (format: TSV)')
parser.add_argument('-cs', '--cases', nargs='+', help='Cases samples IDs')
parser.add_argument('-ct', '--controls', nargs='+', help='Controls samples IDs')
parser.add_argument('-o', '--outlabel', required=True, help='Path label for the output graph')
args = parser.parse_args()

## initialise timer
start = datetime.datetime.now()

## load DataMatrix of cases controls methylation probabilities
inmatrix_df = pd.read_csv(args.inmatrix, sep = "\t")
print("cases controls meth proba df loaded")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## load DataMatrix of testings methylation probabilities
testings_df = pd.read_csv(args.testings, sep = "\t")
print("testings meth proba df loaded")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## get testings spl id
testings_id = [c for c in testings_df.columns.tolist() if c not in ["cpg"] + args.cases + args.controls ]

## store regions in dataframe
reg_df = pd.read_csv(args.sign, sep = "\t", header=None, usecols=[0], names=["cpg"])
print(reg_df)
print("reg dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## restrict cases controls datamatrix on positions from signature
inmatrix_df = inmatrix_df.loc[inmatrix_df['cpg'].isin(reg_df['cpg'].tolist())].reset_index(drop=True)
print(inmatrix_df)
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## restrict testings datamatrix on positions from signature
testings_df = testings_df.loc[testings_df['cpg'].isin(reg_df['cpg'].tolist())].reset_index(drop=True)
print(testings_df)
print("testings_df dataframe restricted")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

##set a common dataframe with methylation probabilities only
full_df = pd.concat([inmatrix_df,testings_df],axis=1)
print(full_df)
print("final dataframe done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)

## prepare plot for correlation
plt.figure(figsize=(6,5), dpi=300)
##plt.title("Pearson correlation of individual PacBio LRS samples\n" + "to the means of Males and Females groups\n", size=12)
plt.plot([0,1] ,[0,1] ,'-', color='black')

## correlate VUS samples
for spl, spl_stat, color in zip(testings_id, ['Female1_test', 'Female2_test', 'Male_test'], ["green", "#d7369e", "gold"]):
    curr_df = full_df
    ## compute mean from cases and control groups
    curr_df['mean_byCpG_controls'] = curr_df.loc[:, args.controls].mean(axis=1)
    curr_df['mean_byCpG_cases'] = curr_df.loc[:, args.cases].mean(axis=1)
    ## compute pearson correlation on testing dataframe
    corr = curr_df[['mean_byCpG_controls', 'mean_byCpG_cases', spl]].corr(method='pearson', numeric_only=True)
    print("Pearson correlation done")
    print(corr)
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
    x = round(corr[spl]['mean_byCpG_cases'], 2)
    y = round(corr[spl]['mean_byCpG_controls'], 2)
    plt.plot(x ,y ,'^', color=color, label=spl_stat)

## correlate on cases and controls while leaving one out
spl_to_iterate = args.cases + args.controls
colors_to_iterate = ["red" for s in args.cases] + ["blue" for s in args.controls]
stat_to_iterate = ["Females_train" for s in args.cases] + ["Males_train" for s in args.controls]
marks_to_iterate = ["s" for s in args.cases] + ["+" for s in args.controls]
for spl_out, color, spl_stat, s_mark in zip(spl_to_iterate, colors_to_iterate, stat_to_iterate, marks_to_iterate) :
    curr_df = full_df
    curr_ctl = [s for s in args.controls if s != spl_out]
    curr_cas = [s for s in args.cases if s != spl_out]
    ## compute mean from cases and control groups
    curr_df['mean_byCpG_controls'] = curr_df.loc[:, curr_ctl].mean(axis=1)
    curr_df['mean_byCpG_cases'] = curr_df.loc[:, curr_cas].mean(axis=1)
    ## compute pearson correlation on testing dataframe
    corr = curr_df[['mean_byCpG_controls', 'mean_byCpG_cases', spl_out]].corr(method='pearson', numeric_only=True)
    print("Pearson correlation done")
    print(corr)
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
    x = round(corr[spl_out]['mean_byCpG_cases'], 2)
    y = round(corr[spl_out]['mean_byCpG_controls'], 2)
    if spl_out == args.cases[0] or spl_out == args.controls[0]:
        plt.plot(x ,y, s_mark, color=color, mfc='none', label=spl_stat)
    else:
        plt.plot(x ,y, s_mark, color=color, mfc='none')
    ##if spl_out in args.cases:
    ##    s_idx = args.cases.index(spl_out)
    ##    plt.annotate("".join([" ", str(s_idx + 1)]), xy=(x, y))

##plt.legend(loc=(1.04, -0.2))
plt.legend(loc=1)
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('Pearson correlation to the mean of Females group', size=12)
plt.ylabel('Pearson correlation to the mean of Males group', size=12)
plt.yticks(size=12)
plt.xticks(size=12)
fig_title = args.outlabel + '.pearson_corr_plot.png'
plt.tight_layout()
plt.savefig(fig_title)
plt.close()
print("Pearson plot done")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)
