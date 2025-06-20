# Ivashchenko-Pacbio-Methylation

Scripts used in analyses from: "_Genome-wide methylation detection and episignature analysis using PacBio long-read sequencing_"

## Downsampling

This script downsamples methylation data for two types:  
1. Illumina EMSeq / MethylSeq (bedGraph format)  
2. PacBio BAM files (long-read data)

**File:** `downsample_methylation_data.sh`  
- For EMSeq: uses `downsample_methylKit.py` (from [methylation_tools](https://github.com/nebiolabs/methylation_tools))  
- For PacBio: uses `samtools view -s` to simulate different read depths (e.g., 5x, 10x, 15x...)  
- Output: downsampled bedGraph and BAM files  

## Differential methylation analysis on PacBio LRS data

Pipeline to identify differentially methylated CpGs and generate a consensus signature.  

**File:** `diff_methylation_pipeline.sh`  
**Steps:**  
1. Format input CpG coverage files into a matrix  
2. Identify differentially methylated CpGs (signature)  
3. Extract consensus CpG signature across iterations  
4. Validate signature using:
   - PCA plot
   - Hierarchical clustering
   - Pearson correlation  

**Input:** `.combined.bed` files for group1, group2, and test samples  
**Output:**  
- Signature CpG lists  
- Train/test matrices  
- Plots and correlation stats  

## Power analysis

Script to model CpG methylation and run power analysis over many sample sizes.

**File:** `power_analysis.sh`  
**Steps:**  
1. Estimate methylation mean and SD using case/control groups  
2. Loop through case/control sizes (10 to 100 by 5)  
3. Run power calculation script on each size combination  
4. Estimate control SD for imprinted regions  
5. Simulate samples (e.g., 500 samples over 20M sites)  
6. Simulate signature detection:
   - Leave one sample out  
   - Mann-Whitney U test  
   - FDR correction  
   - Pearson correlation (observed vs. permuted labels)  

**Output:**  
- Power stats for various sample sizes  
- Simulated methylation data  
- Signature detection accuracy  
