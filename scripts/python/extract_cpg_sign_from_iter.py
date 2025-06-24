#!/usr/bin/env python3

import sys
import argparse
import datetime

###################################
##
## MAIN
##
###################################

## Parameters parser
parser = argparse.ArgumentParser(description="Extract full and consensus CpG signature from Training iterations")
parser.add_argument('-i','--iters', nargs='+', required=True, help='Paths to files with signature CpGs from each iteration (format: TSV)')
parser.add_argument('-o', '--outlabel', required=True, help='Path label for the output files')
args = parser.parse_args()

## initialise timer
start = datetime.datetime.now()

## store cpgs from all iterations
cpg_dict = {}
for file in args.iters:
    print("Parsing " + file)
    print("Time elapsed:", datetime.datetime.now() - start, flush=True)
    with open(file, "r") as it_filin:
        for line in it_filin:
            line = line.strip()
            if line not in cpg_dict:
                cpg_dict[line] = 0
            cpg_dict[line] += 1

## write file with full cpg signature found in at least one iteration
full_sign = args.outlabel + '.full_signature_CpGs.tsv'
## write file with consensus cpg signature common to all iterations
cons_sign = args.outlabel + '.consensus_signature_CpGs.tsv'

with open(full_sign, "w") as filout_full:
    with open(cons_sign, "w") as filout_cons:
        for cg in cpg_dict:
            ##cg_end = int(cg.split("_")[1]) + 1
            ##to_write = cg.split("_") + [str(cg_end)]
            to_write = [cg]
            filout_full.write("\t".join(to_write) + "\n")
            ## write consensus sites
            if cpg_dict[cg] == len(args.iters):
                filout_cons.write("\t".join(to_write) + "\n")
print("Full signature sites " + full_sign)
print("Consensus signature sites " + cons_sign)
print("Time elapsed:", datetime.datetime.now() - start, flush=True)
