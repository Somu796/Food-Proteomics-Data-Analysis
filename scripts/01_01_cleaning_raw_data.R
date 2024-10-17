# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: "data/01_gene_data_post_01_cleaning.RData", 
#           
#             "data/02_sample_details_sample_labels_post_01_cleaning.RData", 
#             

# 1. Library --------------------------------------------------------------
pacman::p_load(tidyverse, janitor, tidylog)

source("scripts/00_path_variables.R")
# 2. Function for data cleaning -------------------------------------------

data_cleaning <- function(rawdata_path, save_gene_data_path, save_sample_details_path){
  # 2. Data Transformation -----------------------------------------------------------------
  
  
  # ST data
  gene_data <- read_csv(rawdata_path)
  
  gene_data <- gene_data %>%
    # remove_col: Cleaning from "Raw abundance" (modify)
    select(all_of(colnames(gene_data)[1:which(colnames(gene_data) == "Raw abundance")-1])) %>%
    # remove_col: From "Peptide count" to "Description", it is not important information (modify)
    select(- colnames(gene_data)[which(gene_data == "Peptide count", arr.ind = TRUE)[2]] : 
             - colnames(gene_data)[which(gene_data == "Description", arr.ind = TRUE)[2]-1]) %>%
    # fill_up: filling first two column upwards (recheck : whether the fill_up operation is not changing other values than the c(1,2) rows)
    fill(c(1,2), .direction = "up") 
  
  gene_data[1,] <- as.data.frame(t(gene_data[1,])) %>% #fill_side: Filling Sample names side-wise
    fill(1, .direction = "down") %>%
    t(.)
  
  # 3. Saving Identifiers (Sample Names and Label Names) ---------------------------------------------------
  
  ## 3.1. Getting the Sample vs Label data frame ------------------------ 
  sample_details <- gene_data[c(2,1),3:ncol(gene_data)] %>% #Keeping Sample and Label Names saved for future
    t() %>%
    as.data.frame()
  
  rownames(sample_details) <- NULL # Converting rownames to NULL
  colnames(sample_details) <- c("sample", "label") # renaming colnames
  
  # 4. Replacing the Column name with 2nd row  ------------------------------
  
  gene_data <- gene_data %>% 
    # rownames_to_colnames: Changing the col names with the second row, because it is unique
    # And column names has to be unique
    row_to_names(row_number = 2)
  
  
  # 5. Saving the data in .RData format for future use ----------------------
  gene_data
  save(gene_data, file = save_gene_data_path)
  sample_details
  save(sample_details, file = save_sample_details_path)
  
}


# 3. Defining data path and data saving places ----------------------------

# save_gene_data_path = c("data/01_gene_data_post_01_cleaning.RData")
# save_sample_details_path = c("data/02_sample_details_sample_labels_post_01_cleaning.RData")

# 4. Running the functions ------------------------------------------------
data_cleaning(rawdata_path, save_gene_data_path, save_sample_details_path)

# 5. Reading generated data -----------------------------------------------
## Run section 3 to individually run it

# load(save_gene_data_path)
# 
# load(save_sample_details_path)

# rm(list = c("gene_data", "sample_details"))



# DISCARDED: If I want to do one by one --------------------------
# # 2. Data Transformation -----------------------------------------------------------------
# source("scripts/00_path_variables.R")
# 
# # ST data
# gene_data <- read_csv("data/raw data/Progenesis_815_Quantifiable_prot_All.csv")
# 
# gene_data <- gene_data %>%
#   # remove_col: Cleaning from "Raw abundance" (modify)
#   select(all_of(colnames(gene_data)[1:which(colnames(gene_data) == "Raw abundance")-1])) %>%
#   # remove_col: From "Peptide count" to "Description", it is not important information (modify)
#   select(- colnames(gene_data)[which(gene_data == "Peptide count", arr.ind = TRUE)[2]] : 
#            - colnames(gene_data)[which(gene_data == "Description", arr.ind = TRUE)[2]-1]) %>%
#   # fill_up: filling first two column upwards (recheck : whether the fill_up operation is not changing other values than the c(1,2) rows)
#   fill(c(1,2), .direction = "up") 
# 
# gene_data[1,] <- as.data.frame(t(gene_data[1,])) %>% #fill_side: Filling Sample names side-wise
#   fill(1, .direction = "down") %>%
#   t(.)
# 
# # 3. Saving Identifiers (Sample Names and Label Names) ---------------------------------------------------
# 
# ## 3.1. Getting the Sample vs Label data frame ------------------------ 
# sample_details <- gene_data[c(2,1),3:ncol(gene_data)] %>% #Keeping Sample and Label Names saved for future
#   t() %>%
#   as.data.frame()
# 
# rownames(sample_details) <- NULL # Converting rownames to NULL
# colnames(sample_details) <- c("sample", "label") # renaming colnames
# 
# # 4. Replacing the Column name with 2nd row  ------------------------------
# 
# gene_data <- gene_data %>% 
#   # rownames_to_colnames: Changing the col names with the second row, because it is unique
#   # And column names has to be unique
#   row_to_names(row_number = 2)
# 
# 
# # 5. Saving the data in .RData format for future use ----------------------
# gene_data
# save(gene_data, file = "data/01_gene_data_post_01_cleaning.RData")
# sample_details
# save(sample_details, file = "data/02_sample_details_sample_labels_post_01_cleaning.RData")
# 
# 
