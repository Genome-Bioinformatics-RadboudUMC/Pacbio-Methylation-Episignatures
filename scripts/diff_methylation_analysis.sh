#!/bin/bash

# Differential Methylation Analysis Script
# Steps:
# 1. Format CpG data into matrix 
# 2+3. Identify CpGs with differential methylation (signature)
# 3. Extract final CpG signature from iteration results
# 4. Validate with PCA, clustering, Pearson correlation

######################
# Step 1: Format Data
######################

format_script=/path/to/scripts/formatCpg_forDataMatrix.py

input_controls=$(ls /path/to/controls/P*.combined.bed)
input_cases=$(ls /path/to/cases/DNA*.combined.bed)
input_new_cases=$(ls /path/to/new_cases/DNA*.combined.bed)

output_label=/path/to/output/DataMatrix_mincov6_allCas_80pcCtls

declare -a m_arr=("Spl_4" "Spl_6" "Spl_1" "Spl_10" "Spl_9" "VUS2" "Spl_12" \
  "Control17" "Control20" "Control22" "Control19" "Control23" "Control24" \
  "Control4" "Control10" "Control2" "Control5" "Control9" "Control8" \
  "Control6" "Control39" "Control27" "Control30" "Control34" "Control36" \
  "Control29" "Control33" "Control37")

declare -a f_arr=("Spl_5" "Spl_3" "Spl_11" "VUS1" \
  "Control13" "Control16" "Control38" "Control18")

declare -a test_arr=("Spl_7" "Control7" "Control14")

python ${format_script} \
  -i ${input_new_cases} ${input_controls} ${input_cases} \
  -cs "${f_arr[@]}" \
  -ct "${m_arr[@]}" \
  -p "${test_arr[@]}" \
  -m 6 -r 0.8 \
  -o ${output_label}

#############################################
# Step 2+3: Run Signature Detection and Output
#############################################

signature_script=/path/to/scripts/getSignature_onDataMatrix_onMedians.py
iter=24

input_matrix=${output_label}_train_data.tsv
metrics_file=${output_label}_cases_controls_metrics_file_onMedians_iter${iter}.tsv

output_diff_matrix=/path/to/output/Cases_Controls_DiffMeansMatrix_iter${iter}.tsv
output_signature=/path/to/output/Cases_Controls_SignatureCpGs_iter${iter}.tsv

python ${signature_script} \
  -i ${input_matrix} \
  -m ${metrics_file} \
  -dmin 20 \
  -cs "${f_arr[@]}" \
  -ct "${m_arr[@]}" \
  -om ${output_diff_matrix} \
  -os ${output_signature} \
  -p 0.05

########################################################
# Step 3: Extract Consensus CpG Signature from Results
########################################################

extract_script=/path/to/scripts/extract_cpg_sign_from_iter.py

signature_files=$(ls /path/to/output/Cases_Controls_SignatureCpGs_*_iter*.tsv)
final_signature_label=/path/to/output/Final_SignatureCpGs

python ${extract_script} \
  -i ${signature_files} \
  -o ${final_signature_label}

##########################
# Step 4: Validation Step
##########################

signature_file=${final_signature_label}.tsv
train_matrix=${output_label}_train_data.tsv
test_matrix=${output_label}_test_data.tsv

# PCA
pca_script=/path/to/scripts/get_PCA_plot.py
pca_out=/path/to/output/PCA_onSignatureCpGs

python ${pca_script} \
  -i ${train_matrix} -ts ${test_matrix} -t ${signature_file} \
  -cs "${f_arr[@]}" -ct "${m_arr[@]}" -tsm "${test_arr[@]}" -o ${pca_out}

# Hierarchical Clustering
clust_script=/path/to/scripts/get_hierarchical_clustering.py
clust_out=/path/to/output/Clustermap_onSignatureCpGs

python ${clust_script} \
  -i ${train_matrix} -ts ${test_matrix} -t ${signature_file} \
  -cs "${f_arr[@]}" -ct "${m_arr[@]}" -tsm "${test_arr[@]}" -g "clustermap" -o ${clust_out}

# Pearson Correlation
pearson_script=/path/to/scripts/pearson_corr_btw_groups.py
pearson_out=/path/to/output/PearsonCorr_15X_8Females_29Males

python ${pearson_script} \
  -i ${train_matrix} -ts ${test_matrix} -s ${signature_file} \
  -cs "${f_arr[@]}" -ct "${m_arr[@]}" -tsm "${test_arr[@]}" -o ${pearson_out}