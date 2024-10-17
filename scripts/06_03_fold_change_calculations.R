# Description
# To RUN: This is an independent script. This has to be run manually. 
# saving_data: yes
# data_name: "data/12_anova_preFold_result_table_post_06.RData"
# This calculates the factorial anova and group by mean for fold change calculation


# To automate: 

# 1. Libraries ------------------------------------------------------------
pacman::p_load(qs, tidyverse, janitor, broom, limma, car, agricolae, glue, tidylog) # rstatix,
source("scripts/00_path_variables.R")


# 2. Calculating the Fold Change ----------------------------------------------
## 2.1. Importing knn_imputed_data -----------------

knn_imputed_data = qs::qread("data_proc.qs") %>% 
  rownames_to_column(var = "sample")

## 2.2. Loading "sample_details_01" and merging with "knn_imputed_data" -------------------

load(sample_details_01) 
knn_imputed_data <- sample_details_01 %>% 
  full_join(knn_imputed_data, by = join_by("sample"))

## 2.3. Getting aggregated mean per group, for "Horn" and "NoHorn" (modify) ------------------

# Calculating aggregated means
fold_change_result = c()

groupby_Stress = aggregate(knn_imputed_data[, (ncol(sample_details_01)+1):ncol(knn_imputed_data)], by=list(knn_imputed_data[[4]]), FUN=mean)
groupby_Stress <- groupby_Stress %>% 
  rename(Group = Group.1)


groupby_Stress_Feeding = aggregate(knn_imputed_data[, (ncol(sample_details_01)+1):ncol(knn_imputed_data)], by=list(knn_imputed_data[[4]], knn_imputed_data[[3]]), FUN=mean)

groupby_Stress_Feeding <- groupby_Stress_Feeding %>% 
  arrange(Group.1) %>% 
  mutate(Group = paste(Group.1, Group.2, sep = "."), .before = 2) %>% 
  select(- c("Group.1", "Group.2"))

fold_change_result = rbind(groupby_Stress, groupby_Stress_Feeding)

# In general, it is a independent fold change
# for (i in 3:4){
#   groupby = aggregate(knn_imputed_data[, (ncol(sample_details_01)+1):ncol(knn_imputed_data)], by=list(knn_imputed_data[[i]]), FUN=mean)
#   fold_change_result = rbind(fold_change_result, groupby)
# }

# Converting them for merging with ANOVA table
fold_change_result = fold_change_result %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  row_to_names(1) %>% 
  mutate_all(as.numeric) %>% 
  rownames_to_column(var = "Accession")

## 5.4. Merging with the TAB0 table to create a new "anova_preFold_result_table" ------------
load(gene_data_02)
fold_change_result <- gene_data_02 %>% 
  select(Accession, gene_name, Protein_name) %>% 
  right_join(fold_change_result, by = "Accession")

save(fold_change_result, file = save_fold_change_result)
