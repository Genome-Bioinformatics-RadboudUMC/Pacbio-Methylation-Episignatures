#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from scipy.stats import mannwhitneyu
import datetime

###################################
##
## MAIN
##
###################################

## Parameters parser
parser = argparse.ArgumentParser(description="Format PacBio CpGs to data matrix and restrain to sites detected in all samples")
parser.add_argument('-i','--inmatrix', required=True, help='Path to cases controls datamatrix with methylation probabilities (format: TSV)')
parser.add_argument('-cs', '--cases', nargs='+', required=True, help='Cases samples IDs')
parser.add_argument('-ct', '--controls', nargs='+', required=True, help='Controls samples IDs')
parser.add_argument('-o', '--outfile', required=True, help='Path to cases controls metrics file')
args = parser.parse_args()

## initialise timer
start = datetime.datetime.now()

## write metrics for cases controls groups to outfile
counter = 0
with open(args.outfile, "w") as filout_metrics:
    header = "cpg" + "\t" + "\t".join(["median_cases", "median_controls", "medians_diff", "SD_cases", "SD_controls", "pval_raw"]) + "\n"
    filout_metrics.write(header)
    print("Created Metrics Output file " + args.outfile)
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
    with open(args.inmatrix, "r") as filin:
        for line in filin:
            line = line.strip().split("\t")
            counter += 1
            if "cpg" in line:
                ## from header of datamatrix get indexing of cases and controls
                cpg_id = line[0]
                cas_idx = {spl: line.index(spl) for spl in args.cases}
                ctl_idx = {spl: line.index(spl) for spl in args.controls}
            else:
                ## count number of lines parsed
                if counter % 100000 == 0:
                    print("100,000 lines parsed")
                    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
                ## iterate on data lines from datamatrix
                cpg = line[0]
                ## get metrics for cases and controls groups
                median_cs = np.median([float(line[cas_idx[s]]) for s in args.cases])
                median_ctl = np.median([float(line[ctl_idx[s]]) for s in args.controls])
                diff_medians = median_cs - median_ctl
                if abs(diff_medians) >= 5:
                    sd_cs = np.std([float(line[cas_idx[s]]) for s in args.cases], ddof=1)
                    sd_ctl = np.std([float(line[ctl_idx[s]]) for s in args.controls], ddof=1)
                    pval = mannwhitneyu([float(line[cas_idx[s]]) for s in args.cases], [float(line[ctl_idx[s]]) for s in args.controls]).pvalue
                    ## append in corresponding iteration file
                    to_write = [cpg] + [str(median_cs), str(median_ctl), str(diff_medians), str(sd_cs), str(sd_ctl), str(pval)]
                    filout_metrics.write("\t".join(to_write) + "\n")
print("Process complete")
print("Time elapsed:", datetime.datetime.now() - start, flush=True)
