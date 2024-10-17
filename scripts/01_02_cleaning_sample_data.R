# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: 
# input
# sample_details
#output
# sample_details_01

# 1. Loading Libraries ----------------------------------------------------
pacman::p_load(MetaboAnalystR, tidyverse, tidylog)
source("scripts/00_path_variables.R")

# 2. Cleaning brackets of sample_details$sample and converting to more descriptive labels --------
load(sample_details)


# 3. Converting to more descriptive labels -----------------------------------
source("scripts/functions/cleanLabels.R") # this to convert "NoStress_NCM -> NoStress"

## 3.1. Define the label mappings for each class ---------------------------
label_mappings <- list(
  "Stress" = list(
    list(pattern = "NoStress", new_label = "NoStress"),
    list(pattern = ".*", new_label = "Stress")
  ),
  "Feeding_Regime" = list(
    list(pattern = "_LEP_", new_label = "Lipid_VitE_PlantExt"),
    list(pattern = "_LC_", new_label = "Lipid_VitE"),
    list(pattern = "_L_", new_label = "Lipid"),
    list(pattern = "_T_", new_label = "Control"),
    list(pattern = ".*", new_label = "NA")
  )
)

## 3.2. Adding the clean labels to the "sample_details" data.frame ---------------------------

sample_details_01 = add_column(sample_details,
                               Feeding_Regime =  as.vector(cleanLabels(sample_details$label, "Feeding_Regime", label_mappings)),
                               Slaughter_Condition =  as.vector(cleanLabels(sample_details$label, "Stress", label_mappings)),
                               .after = 2)

# 4. Saving Details -------------------------------------------------------

save(sample_details_01, file = save_sample_details_01_path)

