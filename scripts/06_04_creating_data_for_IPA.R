# Description
# To RUN: This is an independent script.
# saving_data: YES
# data_name:
# This produces NA in case of post hoc test related proteins

# 1. Loading Libraries ----------------------------------------------------
pacman::p_load(xlsx, tidyverse, glue,tidylog)


# 2. Loading the data -----------------------------------------------------

source("scripts/00_path_variables.R")
source(getFoldChange)

factorial_anova_result <- read.csv(save_factorial_anova_result) # for venn diagram
factorial_anova_posthoc_result <- read.csv(save_factorial_anova_posthoc_result) # for upset plot
load(fold_change_result)
# load(gene_data_02)

# 3. Saving, accession, gene_name, protein name, significance value and fold change ---------------


## 3.1. Slaughter Condition -----------------------

# Data in Anova result
data_for_foldchange_cal <- factorial_anova_result %>%
  full_join(fold_change_result, by = c("Accession", "gene_name","Protein_name")) %>%
  select(- "Residuals")

Factor = "Slaughter_Condition"
meanFactorLevel.1 = "Stress"
meanFactorLevel.2 = "NoStress"

columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
FOldChangeData_Slaughter = c()
FOldChangeData_Slaughter = getFoldChange(data_for_foldchange_cal, columns)

FOldChangeData_Slaughter <- FOldChangeData_Slaughter %>% 
  select(- c(4:6))

colnames(FOldChangeData_Slaughter) [3:4] <- c(paste0("p_",Factor), paste0("log2FC_",Factor))

## 3.2. Feeding Regime -------------------------

# Data in the Post Hoc Result
data_for_foldchange_cal_posthoc <- factorial_anova_posthoc_result %>% 
  full_join(fold_change_result, by = c("Accession", "gene_name","Protein_name"))

###  Column name processing -------

#### Separating the No Stress and Stress df separately -------------------------------
nostress_colnames <- grepl("NoStress", colnames(data_for_foldchange_cal_posthoc))
stress_colnames <- setdiff(colnames(data_for_foldchange_cal_posthoc), colnames(data_for_foldchange_cal_posthoc)[nostress_colnames])
nostress_colnames[1] <- TRUE # Accession
nostress_colnames[2] <- TRUE # gene_name

#### Modifying the column names -----------------------

##### replacing "..." with "-"  ----------------------
# nostress_df
nostress_df <- data_for_foldchange_cal_posthoc[,nostress_colnames]
colnames(nostress_df) <- str_replace(colnames(nostress_df),"^NoStress\\.\\.\\.", "") %>% 
  str_replace(., "\\.\\.\\."," \\- ")

colnames(nostress_df) <- str_replace(colnames(nostress_df),"^NoStress\\.", "")

# stress_df
stress_df <- data_for_foldchange_cal_posthoc[,stress_colnames]
colnames(stress_df) <- str_replace(colnames(stress_df),"^Stress\\.\\.\\.", "") %>% 
  str_replace(., "\\.\\.\\."," \\- ")
colnames(stress_df) <- str_replace(colnames(stress_df),"^Stress\\.", "")

#####  swapping "stressA-stressB" to "stressB-stressA" ---------------------------

swap_words <- function(input_vector) {
  new_vector <- c()
  for(i in 1: length(input_vector)){
    modified_input_string <- str_replace(input_vector[i], "(\\w+) - (\\w+)", "\\2 - \\1")
    new_vector <- append(new_vector, modified_input_string)
  }
  return(new_vector)
}

colnames(nostress_df) <- swap_words(colnames(nostress_df))
colnames(stress_df) <- swap_words(colnames(stress_df))


### 3.2.1. nostress_df ------------------
colnames(nostress_df)

Factor = c("Lipid_VitE_PlantExt - Lipid_VitE", 
           "Lipid_VitE_PlantExt - Lipid",
           "Lipid_VitE - Lipid")

meanFactorLevel.1 = c("Lipid_VitE_PlantExt", 
                      "Lipid_VitE_PlantExt",
                      "Lipid_VitE")

meanFactorLevel.2 = c("Lipid_VitE", 
                      "Lipid",
                      "Lipid")


# Factor = c("Lipid_VitE_PlantExt - Lipid_VitE", 
#            # "Lipid - Control", 
#            # "Lipid_VitE_PlantExt - Control", 
#            "Lipid_VitE_PlantExt - Lipid",
#            # "Lipid_VitE - Control"
#            )
# 
# meanFactorLevel.1 = c("Lipid_VitE_PlantExt", 
#                       # "Lipid",
#                       # "Lipid_VitE_PlantExt",
#                       "Lipid_VitE_PlantExt",
#                       # "Lipid_VitE"
#                       )
# 
# meanFactorLevel.2 = c("Lipid_VitE", 
#                       # "Control",
#                       # "Control",
#                       "Lipid",
#                       # "Control"
#                       )

nostress_IPA <- nostress_df[,1:2]

for (i in 1:length(Factor)){
  
  columns = c("Accession", "gene_name", Factor[i], meanFactorLevel.1[i], meanFactorLevel.2[i])
  FOldChangeData_Fed <- c()
  FOldChangeData_Fed <- nostress_df %>%
    dplyr::select(all_of(columns)) %>%
    filter(., !is.na(.[[Factor[i]]])) 
  
  FOldChangeData_Fed$MeanRatio = FOldChangeData_Fed[[meanFactorLevel.1[i]]]/FOldChangeData_Fed[[meanFactorLevel.2[i]]] #Calculating Mean ratio
  FOldChangeData_Fed$log2FC = log2(FOldChangeData_Fed$MeanRatio)
  
  FOldChangeData_Fed <- FOldChangeData_Fed%>% 
    select(- c(4:6))
  
  colnames(FOldChangeData_Fed) [3:4] <- c(paste0("p_",Factor[i]), paste0("log2FC_",Factor[i]))
  
  nostress_IPA <- nostress_IPA %>% 
    full_join(FOldChangeData_Fed, by = c("Accession", "gene_name"))
}

colnames(nostress_IPA)[-1:-2] <- paste0("nostress_", colnames(nostress_IPA)[-1:-2])

### 3.2.1. stress_df ------------------

colnames(stress_df)

Factor = c("Lipid_VitE_PlantExt - Lipid_VitE", 
           "Lipid_VitE_PlantExt - Lipid",
           "Lipid_VitE - Lipid")

meanFactorLevel.1 = c("Lipid_VitE_PlantExt", 
                      "Lipid_VitE_PlantExt",
                      "Lipid_VitE")

meanFactorLevel.2 = c("Lipid_VitE", 
                      "Lipid",
                      "Lipid")

stress_IPA <- stress_df[,1:2]

for (i in 1:length(Factor)){
  
  columns = c("Accession", "gene_name", Factor[i], meanFactorLevel.1[i], meanFactorLevel.2[i])
  FOldChangeData_Fed <- c()
  FOldChangeData_Fed <- stress_df %>%
    dplyr::select(all_of(columns)) %>%
    filter(., !is.na(.[[Factor[i]]])) 
  
  FOldChangeData_Fed$MeanRatio = FOldChangeData_Fed[[meanFactorLevel.1[i]]]/FOldChangeData_Fed[[meanFactorLevel.2[i]]] #Calculating Mean ratio
  FOldChangeData_Fed$log2FC = log2(FOldChangeData_Fed$MeanRatio)
  
  FOldChangeData_Fed <- FOldChangeData_Fed%>% 
    select(- c(4:6))
  
  colnames(FOldChangeData_Fed) [3:4] <- c(paste0("p_",Factor[i]), paste0("log2FC_",Factor[i]))
  
  stress_IPA <- stress_IPA %>% 
    full_join(FOldChangeData_Fed, by = c("Accession", "gene_name"))
}


colnames(stress_IPA)[-1:-2] <- paste0("stress_", colnames(stress_IPA)[-1:-2])

# nostress_IPA
# stress_IPA
# FOldChangeData_Slaughter

IPA_Data <- FOldChangeData_Slaughter %>% 
  full_join(nostress_IPA, by = c("Accession", "gene_name")) %>% 
  full_join(stress_IPA, by = c("Accession", "gene_name"))

# 4. Saving data for IPA Analysis -----------------------------------------

write.csv(IPA_Data, file = IPA_Data_dest, row.names = FALSE)

# ####  Lipid_VitE_PlantExt - Lipid_VitE
# 
# Factor = "Lipid_VitE_PlantExt - Lipid_VitE"
# meanFactorLevel.1 = "Lipid_VitE_PlantExt"
# meanFactorLevel.2 = "Lipid_VitE"
# 
# columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
# FOldChangeData_Fed_1 = c()
# FOldChangeData_Fed_1 = getFoldChange(nostress_df, columns) 
# FOldChangeData_Fed_1 <- FOldChangeData_Fed_1%>% 
#   select(- c(4:6))
# colnames(FOldChangeData_Fed_1) [3:4] <- c(paste0("p_",Factor), paste0("log2FC_",Factor))
# 
# 
# 
# #### Lipid - Control 
# 
# Factor = "Lipid - Control"
# meanFactorLevel.1 = "Lipid"
# meanFactorLevel.2 = "Control"
# 
# columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
# FOldChangeData_Fed_2 = c()
# FOldChangeData_Fed_2 = getFoldChange(nostress_df, columns)
# FOldChangeData_Fed_2 <- FOldChangeData_Fed_2%>% 
#   select(- c(4:6))
# 
# colnames(FOldChangeData_Fed_2) [3:4] <- c(paste0("p_",Factor), paste0("log2FC_",Factor))
# 
# #### Lipid_VitE_PlantExt - Control 
# 
# Factor = "Lipid_VitE_PlantExt - Control"
# meanFactorLevel.1 = "Lipid_VitE_PlantExt"
# meanFactorLevel.2 = "Control"
# 
# columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
# FOldChangeData_Fed_3 = c()
# FOldChangeData_Fed_3 = getFoldChange(nostress_df, columns)
# FOldChangeData_Fed_3 <- FOldChangeData_Fed_3 %>% 
#   select(- c(4:6))
# 
# colnames(FOldChangeData_Fed_3) [3:4] <- c(paste0("p_",Factor), paste0("log2FC_",Factor))
# 
# #### Lipid_VitE_PlantExt - Lipid
# 
# Factor = "Lipid_VitE_PlantExt - Lipid"
# meanFactorLevel.1 = "Lipid_VitE_PlantExt"
# meanFactorLevel.2 = "Lipid"
# 
# columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
# FOldChangeData_Fed_4 = c()
# FOldChangeData_Fed_4 = getFoldChange(nostress_df, columns)
# FOldChangeData_Fed_4 <- FOldChangeData_Fed_4 %>% 
#   select(- c(4:6))
# 
# colnames(FOldChangeData_Fed_4) [3:4] <- c(paste0("p_",Factor), paste0("log2FC_",Factor))
# 
# #### Lipid_VitE - Control
# 
# Factor = "Lipid_VitE - Control"
# meanFactorLevel.1 = "Lipid_VitE"
# meanFactorLevel.2 = "Control"
# 
# columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
# FOldChangeData_Fed_5 = c()
# FOldChangeData_Fed_5 = getFoldChange(nostress_df, columns)
# FOldChangeData_Fed_5 <- FOldChangeData_Fed_5 %>% 
#   select(- c(4:6))
# 
# colnames(FOldChangeData_Fed_5) [3:4] <- c(paste0("p_",Factor), paste0("log2FC_",Factor))
# 
# nostress_IPA <-  FOldChangeData_Fed_1 %>% 
#   full_join(FOldChangeData_Fed_2, by = c("Accession", "gene_name")) %>% 
#   full_join(FOldChangeData_Fed_3, by = c("Accession", "gene_name")) %>% 
#   full_join(FOldChangeData_Fed_4, by = c("Accession", "gene_name")) %>% 
#   full_join(FOldChangeData_Fed_5, by = c("Accession", "gene_name"))
# 




