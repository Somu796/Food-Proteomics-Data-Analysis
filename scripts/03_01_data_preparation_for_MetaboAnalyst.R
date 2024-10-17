# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: 
# "data/08_gene_data_post_removing_NA_entries_from_gene_primary_after_03_01_script.RData"
# "data/09_data_for_metabo_analysis_normalization_after_03_01.RData"
# "data/09_data_for_metabo_analysis_normalization_after_03_01.csv"

# Very Important: Here, we are filtering out all the "Protein_name" entries which has corresponding "Gene_primary" as NA. 
# Also after modifying the duplicated gene names, the data frame is being saved as gene_data_02

# 1. Libraries -----------------------------------------------------------

pacman::p_load(tidyverse, janitor, naniar, tidylog)
source("scripts/00_path_variables.R")

## Input Data
# enriched_gene_data_01 <- c("data/07_enriched_gene_data_01_after_02_02_merging_gene_data_with_uniprot.RData")
# sample_details_path <- c("data/02_sample_details_sample_labels_post_01_cleaning.RData")

## Output Data
# save_gene_data_02_path <- c("data/08_gene_data_post_removing_NA_entries_from_gene_primary_after_03_01_script.RData")
# save_data_for_metabo_path <- c("data/09_data_for_metabo_analysis_normalization_after_03_01.RData")
# save_csv_data_for_metabo_path <- c("data/09_data_for_metabo_analysis_normalization_after_03_01.csv")

# 2. Running in a loop ----------------------------------------------------

# 2. Load "enriched_gene_data_01" and modify duplicated gene names to make it unique ------------------

load(enriched_gene_data_01) # enriched gene data
gene_data_02 <- enriched_gene_data_01[!is.na(enriched_gene_data_01$Gene_primary),] #removing all the entries having "Gene_primary" as NA, means only protein names were not considered

dup_indices <- c()
dup_indices <- duplicated(gene_data_02$Gene_primary) | duplicated(gene_data_02$Gene_primary, fromLast = TRUE) # index of duplicate values
# the above OR logic has been used to tackle the limitation on duplicated function, 
# eg. a = c(6, 6, 8, 8) \n duplicated(a) -> F, T, F, T \n duplicated(a, fromLast = TRUE) ->  T, F, T ,F; 
# to fix this OR logic will give back all the repetition, ask chatgpt

# Add .1, .2, etc., to repeated values
gene_data_02$gene_name <- gene_data_02$Gene_primary
gene_data_02$gene_name[dup_indices] <- ave(as.vector(gene_data_02$Gene_primary[dup_indices]), 
                                           as.vector(gene_data_02$Gene_primary[dup_indices]), 
                                           FUN = function(x) paste0(x, ".", seq_along(x)))

gene_data_02 <- gene_data_02 %>% #relocating the new created column "gene_name", which has been created to handle the duplicate gene names 
  #in "Gene_primary", DMD, DMD -> DMD.1, DMD.2  
  relocate(gene_name, .after = 1)


# 3. Modifying data for MetaboAnalyst (here unnecessary columns for data pre-processing were removed) -------------------------------------
data_for_metabo <-  gene_data_02 %>% 
  #only keeping "gene_name"
  select(c(-2:-8)) %>%  # select(c(-1, -3:-8)) for ***gene_name*** as col name----*************************************
  #transporting gene names to columns
  t() %>% 
  as.data.frame() %>%
  # first row is made as the col names
  row_to_names(1) %>% 
  # row_names to first column
  rownames_to_column("sample")


# 4. Merging "data_for_metabo" with sample labels to get the Labels -----------
load(sample_details)
data_for_metabo <- sample_details %>% 
  # adding labels with the "data_for_metabo"
  full_join(data_for_metabo, by = join_by("sample")) %>% 
  rename(all_of(c(Sample = "sample", Label = "label")))

data_for_metabo[data_for_metabo==0] <- NA #converting all zeros to NA, for easy counting 


# 5. Missing data (%) in the "data_for_metabo" ---------------------------------------------

## 5.1. Misising Data table ----------------
missing_data_perc <- 
  # This gives the missing value percentage
  colMeans(is.na(data_for_metabo)) %>%
  # Converting to dataframe
  as.data.frame() %>% 
  # "gene_name" are in rowname, getting in the column
  rownames_to_column(var = "gene_name") %>% 
  # While converting to data.frame "missing_perc" came as ".", converting that
  rename(., missing_perc = `.`) %>% 
  # Ordering missing data in the decreasing order
  arrange(desc(missing_perc)) %>%
  # Filtering out the gene_name with 0 missing data
  filter(missing_perc != 0) %>% 
  # transposing it 
  t() %>% 
  # converting to data frame
  as.data.frame() %>%
  # gene names are converted as column names
  row_to_names(1)


## 5.2. Missing Data Over All Hetamap -----
df_sorted <- data_for_metabo %>%
  # select(-"Label") %>% 
  select(order(colSums(is.na(.)), decreasing = TRUE))

missing_cols <- colnames(df_sorted)[colSums(is.na(df_sorted)) > 0]

missing_cols_counts <- colSums(is.na(df_sorted))
nonzero_missing_counts <- missing_cols_counts[missing_cols_counts > 0] 

max_missing_index <- which.max(missing_cols_counts) # to find index of max value
min_missing_index <- which(names(missing_cols_counts) == names(which.min(rev(nonzero_missing_counts)))) # to find index of min values
# all missing value columns == column which have minimum, after reversing column names because if two columns have the min, it will only consider the first not last

xmin <- max_missing_index - 0.5
xmax <- min_missing_index + 0.5
ymin <- 0.5
ymax <- nrow(df_sorted) + 0.5

missing_data_hetamap_overall_01 <- vis_miss(df_sorted) + 
  labs(x = "Genes") +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "red", linewidth = 1) + 
  labs(x = "Genes") +
  theme(axis.text.x = element_text(size = 0.1))+
  scale_y_continuous(breaks = seq(1, length(rownames(df_sorted)), 1))

# print(missing_data_hetamap_overall_01)


## 5.3. Missing Data Specific ----- 

df_sorted_specific <-  df_sorted %>% 
  select(all_of(missing_cols))

missing_data_hetamap_specfic_02 <- vis_miss(df_sorted_specific, show_perc = FALSE) + 
  labs(x = "Genes") +
  theme(axis.text.x = element_text(size = 8, face = "bold", angle=90))+
  scale_y_continuous(breaks = seq(1, length(rownames(df_sorted_specific)), 1))

missing_data_hetamap_specfic_02


# 6. Saving Data -----------------------------------------------------------

save(gene_data_02, file = save_gene_data_02_path)
save(data_for_metabo, file = save_data_for_metabo_path)
write.csv(data_for_metabo, file = save_csv_data_for_metabo_path, row.names = FALSE)


missing_data_perc
missing_data_hetamap_overall_01
missing_data_hetamap_specfic_02

