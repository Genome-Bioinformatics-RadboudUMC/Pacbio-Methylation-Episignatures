#!/bin/bash

# Step 1: Model CpG methylation levels
# - Calculates mean and SD for case and control groups
# - Used later to simulate realistic CpG data
Rscript modelCpGMethylation.R

# Step 2: Run power calculations across a range of sample sizes
# - Loops through combinations of case and control sample counts
# - Runs power analysis for each combo
for case_samples in {10..100..5}; do
  for control_samples in {10..100..5}; do
    echo "Running calculatePower.R with case_samples=$case_samples, control_samples=$control_samples"
    Rscript scripts/R/calculatePower.R $case_samples $control_samples
  done
done

# Step 3: Estimate minimum standard deviation for PacBio controls
# - Uses imprinted regions (bed file)
Rscript pacbioControlSDMins.R data/imprinted_regions.txt $samples 100

# Step 4: Simulate samples
# - Simulates 500 samples across 20 million sites
Rscript simulateSamples.R 500 20000000

# Step 5: Simulate methylation signatures
# - Leaves out one sample for testing
# - Runs Mann-Whitney U test
# - Applies FDR correction
# - Computes observed vs. expected Pearson correlations
Rscript simulateSignatures.R