# Ivashchenko-Pacbio-Methylation

Repository with scripts used in analyses from: "_Genome-wide methylation detection and episignature analysis using PacBio long-read sequencing_"

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
2. Identify differentially methylated CpGs (signature) using median-based stats  
3. Extract consensus CpG signature across iterations  
4. Validate signature using:
   - PCA plot
   - Hierarchical clustering
   - Pearson correlation  

**Input:** `.combined.bed` files for controls, cases, and test samples  
**Output:**  
- Signature CpG lists  
- Train/test matrices  
- Plots and correlation stats  

## Power analysis
