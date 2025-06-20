#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import datetime

###################################
##
## MAIN
##
###################################

## Parameters parser
parser = argparse.ArgumentParser(description="Format PacBio CpGs to data matrix and restrain to sites detected in all samples")
parser.add_argument('-i','--inputs', nargs='+', required=True, help='Path to PacBio methylation combined detection files (format: BED)')
parser.add_argument('-h1','--hap1', nargs='+', help='Path to PacBio methylation hap1 detection files (format: BED)')
parser.add_argument('-h2','--hap2', nargs='+', help='Path to PacBio methylation hap2 detection files (format: BED)')
parser.add_argument('-cs', '--cases', nargs='+', required=True, help='Cases samples IDs')
parser.add_argument('-ct', '--controls', nargs='+', required=True, help='Controls samples IDs')
parser.add_argument('-p', '--predictset', nargs='+', required=True, help='Samples IDs with not known status')
parser.add_argument('-m', '--mincov', type=int, help='Set a minimal coverage threshold to filter cpg sites detection')
parser.add_argument('-r', '--ctlratio', type=float, help='Minimal percentage of controls that have to fit the mincov threshold')
parser.add_argument('-o', '--outlabel', required=True, help='Path to path label for train and test output data matrix files')
args = parser.parse_args()

## initialise timer
start = datetime.datetime.now()

## restrict PacBio methylation detection to cpg in chr_list and detected in all samples
chr_list = ["chr" + str(nb) for nb in range(1,23)] + ["chrX"]
cgPos_byfile = {}
cov_onCtl = {}
spl_l = []
for bed_f in args.inputs:
    spl = bed_f.split("/")[-1].split(".")[0]
    selected_spls = args.cases + args.controls + args.predictset
    if spl in selected_spls:
        spl_l.append(spl)
        print("Parsing file " + bed_f)
        print("Time elapsed:", datetime.datetime.now() - start, flush=True)
        with open(bed_f, "r") as filin:
            for line in filin:
                line = line.strip().split("\t")
                chr = line[0]
                cpg_start = line[1]
                cpg_pb = line[3]
                cov = int(line[5])
                cov_check = "PASS"
                if chr in chr_list:
                    curr_cgID = chr + "_" + cpg_start
                    ## control coverage on cases
                    if args.mincov is not None:
                        if spl in args.cases:
                            cov_check = "PASS" if cov >= args.mincov else "SKIP"
                        elif spl in args.controls:
                            ## store coverage values on controls
                            if curr_cgID not in cov_onCtl:
                                cov_onCtl[curr_cgID] = []
                            cov_onCtl[curr_cgID].append(cov)
                    if cov_check == "PASS":
                        if curr_cgID not in cgPos_byfile :
                            cgPos_byfile[curr_cgID] = {}
                        cgPos_byfile[curr_cgID][spl] = cpg_pb

## write separately training set and testing/prediction set
train_data = args.outlabel + "_train_data.tsv"
test_data = args.outlabel + "_test_data.tsv"
trainset = args.cases + args.controls
with open(train_data, "w") as filout_train:
    header = "cpg" + "\t" + "\t".join(trainset) + "\n"
    filout_train.write(header)
    with open(test_data, "w") as filout_test:
        header = "cpg" + "\t" + "\t".join(args.predictset) + "\n"
        filout_test.write(header)
        for cpg in cgPos_byfile:
            ## restrict to sites reported by all samples
            if len(cgPos_byfile[cpg]) == len(spl_l):
                ## restrict to sites reported with expected coverage in specified percentage of controls
                if sum([1 for s in cov_onCtl[cpg] if s >= args.mincov])/len(args.controls) >= args.ctlratio:
                    to_write = [cpg] + [cgPos_byfile[cpg][s] for s in trainset]
                    filout_train.write("\t".join(to_write) + "\n")
                    to_write = [cpg] + [cgPos_byfile[cpg][s] for s in args.predictset]
                    filout_test.write("\t".join(to_write) + "\n")
print("Train Data Output Done " + train_data)
print("Test Data Output Done " + test_data)
print("Time elapsed:", datetime.datetime.now() - start, flush=True)
