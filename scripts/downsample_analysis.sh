#!/bin/bash

# This script downsamples methylation data from:
# 1. Illumina MethylSeq/EMSeq (bedGraph format)
# 2. PacBio BAM files (long-read data)

###########################
# Part 1: MethylSeq/EMSeq #
###########################

# Tool: downsample_methylKit.py (from https://github.com/nebiolabs/methylation_tools)
# Goal: Reduce EMSeq bedGraph to ~30x coverage

script=downsample_methylKit.py
input_bg="epiqc_data/EMSeq_HG002_LAB01_LAB02.combined.bedGraph"
output_bg="epiqc_data/EMSeq_HG002_LAB01_LAB02.combined.downsampled30X.bedGraph"

# 0.57 is the fraction estimated to give ~30x
grep -v track "$input_bg" | python "$script" --fraction 0.57 --bedGraph > "$output_bg"

#####################################
# Part 2: PacBio BAM Downsampling   #
#####################################

# Input BAM and output settings
input_bam="/path/to/sample.bam"
output_dir="/path/to/output"
threads=4

# Coverage depths to downsample to
depths=(5 10 15 20 25 30)

# Estimate current average depth
current_coverage=$(samtools depth -a -@ "$threads" "$input_bam" | awk '{sum+=$3} END {print sum/NR}')

# Downsample loop
for target_coverage in "${depths[@]}"; do
  # Skip if current depth is already lower than target
  if (( $(echo "$current_coverage < $target_coverage" | bc -l) )); then
    continue
  fi

  fraction=$(echo "$target_coverage / $current_coverage" | bc -l)
  output_bam="$output_dir/sample.downsampled_${target_coverage}x.bam"

  samtools view -s "$fraction" -@ "$threads" -b "$input_bam" > "$output_bam"
  samtools index -@ "$threads" "$output_bam"
done