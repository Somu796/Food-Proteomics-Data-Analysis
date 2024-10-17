##
# Details: save_gene_data for saving the data frame at the destination, gene_data is for creating a data frame and calling it from destination
# Saving and creating variable can't have same name because, 
# gene_data <- "destination"
# gene_data <- c(1,2,3,4)
# save(gene_data, file = gene_data)
# This doesn't make sense, one name can't direct to two values

# Raw Data
rawdata_path <- "data/raw data/Progenesis_LT_817_Quantifiable_prot_All.csv"

# 2. Datasets -------------------------------------------------------------
## 01 --------
gene_data <- save_gene_data_path <- "data/01_01_gene_data_post_01_cleaning.RData"
sample_details <- save_sample_details_path <- "data/01_02_sample_details_sample_labels_post_01_01_cleaning.RData"
sample_details_01 <- save_sample_details_01_path <- "data/01_03_sample_details_sample_labels_post_01_02_cleaning.RData"
## 02 ------------
### 02_01 -------------
UniprotNames <- save_UniprotNames_path <- "data/02_01_01_enriched_gene_name_information_from_uniprot_post_02_01_enriching_.RData"
merged_Accession_UniprotNames_information_path <- save_merged_Accession_UniprotNames_information_path <- "data/02_01_02_merged_accession_enriched_gene_name_information.RData"
gene_data_01 <- save_gene_data_01_path <- "data/02_01_03_gene_data_post_02_01_with_multple_accession_entries.RData"

### 02_02 -------------
enriched_gene_data_01 <- save_enriched_gene_data_01_path <- "data/02_02_01_enriched_gene_data_01_after_02_02_merging_gene_data_with_uniprot.RData"
unused_accession <- save_unused_accession_path <- "data/02_02_02_unused_accessions_for_enriching_gene_data.RData"

### 03_01 ---------------
gene_data_02 <- save_gene_data_02_path <- "data/03_01_01_gene_data_post_removing_NA_entries_from_gene_primary_after_03_01_script.RData"
data_for_metabo <- save_data_for_metabo_path <- "data/03_01_02_data_for_metabo_analysis_normalization_after_03_01.RData"
csv_data_for_metabo <- save_csv_data_for_metabo_path <- "data/03_01_03_data_for_metabo_analysis_normalization_after_03_01.csv"

### 03_02 ------------------
#Always check "sample" column got modified or not. e.g., remember metaboAnalyst modify"(" <- ""
metabo_data <- save_metabo_data <- "data/03_02_01_metabo_data_for_statistical_analysis_and_enrichment_post_03_02.RData" 

### 03_03 -------------------
residual_model_feeding_regime_knn_imputed <- save_residual_model_feeding_regime_knn_imputed <-  "data/03_03_01_residual_model_feeding_regime_on_knn_imputed_data.csv"
residual_model_feeding_regime_metabo <- save_residual_model_feeding_regime_metabo <- "data/03_03_01_residual_model_feeding_regime_on_metabo_data.csv"
  
### 05_01 ------------------
residual_model_feeding_regime <- save_residual_model_feeding_regime <-  "data/05_01_01_residual_model_feeding_regime.csv"
efs_result_dest <- "data/05_01_02_Ensemble_Feature_selection_output.csv"

### 05_02 ------------------
save_heatmap_EFS_features_Slaughter_Stress <- "plots/05_02_heatmap_EFS_features_Slaughter_Stress.png"

### 06_01 ---------------------

# This has to be added

### 06_02 ------------------
# in case of csv data we don't need factorial_anova_result, only save_.... is needed.
factorial_anova_result <- save_factorial_anova_result <- "data/06_02_01_factorial_anova_result.csv"
factorial_anova_posthoc_result <- save_factorial_anova_posthoc_result <- "data/06_02_02_factorial_anova_posthoc_result.csv"

## 06_03 --------------------
fold_change_result <- save_fold_change_result <- "data/06_03_01_fold_change_result.RData"

## 06_04 --------------------
IPA_Data_dest <- "data/06_04_01_IPA_data_for_analysis.csv"

## 07_01 ----------------------
save_venn_diagram_interaction <- "plots/07_01_01_venn_diagram_interaction.png"
save_venn_diagram_wo_interaction <- "plots/07_01_01_enn_diagram_wo_interaction.png"

save_upset_plot_diet_givenstress <- "plots/07_01_02_upset_plot_diet_givenstress.png"
save_upset_plot_diet_givennostress <- "plots/07_01_02_upset_plot_diet_givennostress.png"

save_volcanoplot_slaughter_condition <- "plots/07_01_03_volcanoplot_slaughter_condition.png"

save_up_and_down_regulated_genes <- "data/07_01_01_up_and_down_regulated_genes.xlsx" 

## 07_02 ----------------------
save_volcano_plot_EFS <- "plots/07_02_01_volcano_plot_EFS.png"
save_EFS_vs_ANOVA <- "plots/07_02_01_volcano_plot_EFS_vs_ANOVA.png"

## 07_03 ------------------------
save_heatmap_significant_genes_Slaughter_Stress_removing_affect_of_factor <- "plots/07_03_01_heatmap_significant_genes_Slaughter_Stress_removing_affect_of_factor.png" # Not This
save_heatmap_significant_genes_Slaughter_Stress_with_metabo_data <- "plots/07_03_02_heatmap_significant_genes_Slaughter_Stress_with_metabo_data.png" # Go with this
save_heatmap_up_and_down_regulated_genes_Slaughter_Stress_with_metabo_data <- "plots/07_03_03_heatmap_up_and_down_regulated_genes_Slaughter_Stress_with_metabo_data.png"
save_heatmap_significant_genes_Feeding_Regime_with_metabo_data <- "plots/07_03_04_heatmap_significant_genes_Feeding_Regime_with_metabo_data.png"
# 3. Functions ----------------------------------------------------------

getAccession <- "scripts/functions/getAccession.R" # 02_01
getGeneNames <- "scripts/functions/getGeneNames.R" # 02_01
gerMergedGeneNames <- "scripts/functions/gerMergedGeneNames.R" #02_01
countAccession <- "scripts/functions/countAccession.R" #02_01
getGeneNamesfromMultipleAccession <- "scripts/functions/getGeneNamesfromMultipleAccession.R" #02_01

getFoldChange = "scripts/functions/getFoldChange.R" #07
getVolcanoPlot = "scripts/functions/getVolcanoPlot.R" #07
changingVennFont = "scripts/functions/changingVennFont.R" #07
