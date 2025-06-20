#!/bin/bash

# This script downsamples methylation data for:
# 1. Illumina MethylSeq/EMSeq (bedGraph format)
# 2. PacBio BAM files (long-read data)

###########################
# Part 1: MethylSeq/EMSeq #
###########################

# Source: https://github.com/nebiolabs/methylation_tools
# Tool: downsample_methylKit.py (must be in PATH or same folder)
# Goal: Downsample EMSeq bedGraph file to approx. 30x

SCRIPT=downsample_methylKit.py
INPUT_BG="epiqc_data/EMSeq_HG002_LAB01_LAB02.combined.bedGraph"
OUTPUT_BG="epiqc_data/EMSeq_HG002_LAB01_LAB02.combined.downsampled30X.bedGraph"

# Fraction (0.57) estimated to reach ~30x coverage
grep -v track "$INPUT_BG" | python "$SCRIPT" --fraction 0.57 --bedGraph > "$OUTPUT_BG"

################################
# Part 2: PacBio BAM Downsampling #
################################

# Inputs
INPUT_BAM="/path/to/sample.bam"
OUTPUT_DIR="/path/to/output"
THREADS=4
DEPTHS=(5 10 15 20 25 30)

# Calculate current average coverage
CURRENT_COVERAGE=$(samtools depth -a -@ "$THREADS" "$INPUT_BAM" | awk '{sum+=$3} END {print sum/NR}')

# Downsample to target depths
for TARGET_COVERAGE in "${DEPTHS[@]}"; do
  if (( $(echo "$CURRENT_COVERAGE < $TARGET_COVERAGE" | bc -l) )); then
    continue
  fi

  FRACTION=$(echo "$TARGET_COVERAGE / $CURRENT_COVERAGE" | bc -l)
  OUTPUT_BAM="$OUTPUT_DIR/sample.downsampled_${TARGET_COVERAGE}x.bam"

  samtools view -s "$FRACTION" -@ "$THREADS" -b "$INPUT_BAM" > "$OUTPUT_BAM"
  samtools index -@ "$THREADS" "$OUTPUT_BAM"
done