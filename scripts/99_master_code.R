
# Master Script
# Only change the location of csv file from  00_path_variables.R

source("scripts/01_01_cleaning_raw_data.R")

rm(list = ls())
source("scripts/01_02_cleaning_sample_data.R")

rm(list = ls())
source("scripts/02_01_enriching_data_with_uniprot_and_description.R")

rm(list = ls())
source("scripts/02_02_merging_gene_data_with_uniprot_data_03.R")

rm(list = ls())
source("scripts/03_01_data_preparation_for_MetaboAnalyst.R")

rm(list = ls())
#Always check "sample" column got modified or not. e.g., remember metaboAnalyst modify"(" <- ""
source("scripts/03_02_data_normalization_with_MetaboAnalyst.R")

rm(list = ls())
# Specific to lipivmus removing variation due to one of the factor
source("scripts/03_03_optional_script_removing_variation_due to_a_factor.R")

rm(list = ls())
# Done only with knn imputed data, variation due to factor Feeding_Regime was not removed
source("scripts/04_01_unsupervised_learning_on_03_02_normalized_data.R")

rm(list = ls())
source("scripts/05_01_supervised_learning_biomarker_prediction.R") # changed things manually for lipivimus

rm(list = ls())
source("scripts/05_02_supervised_learning_biomarker_prediction_plots.R", echo = TRUE)

rm(list = ls())
source("scripts/06_01_assumption_checking_for_ANOVA.R")

rm(list = ls())
source("scripts/06_02_factorial_anova_with_norm_data.R")

rm(list = ls())
source("scripts/06_03_fold_change_calculations.R")

rm(list = ls())
source("scripts/06_04_creating_data_for_IPA.R")

rm(list = ls())
source("scripts/07_01_venndiagram_upSet_volcano_plot.R", echo = TRUE)


rm(list = ls())
source("scripts/07_02_supervised_learning_vs_anova_plots.R", echo = TRUE)


# Only problem with this is the projects has to be saved to display in quarto file,
# but what if the master script is a quarto file,
# plot everything then remove the previous object ---- try it