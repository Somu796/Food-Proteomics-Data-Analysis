# Manual

# 1. install and import the package ---------------------------------------

# install.packages("EFS")
pacman::p_load(tidyverse, glue, janitor , EFS, tidylog)
source("scripts/00_path_variables.R")

# 2. Importing Data ----------------------------------------------

# sample details
load(sample_details_01)

## 2.1. for knn_Imputed Data --------------------------
# sample data
knn_imputed_data <- qs::qread("data_proc.qs")%>% 
  rownames_to_column(var = "sample")

# Data merging
Data_knn_imputed <- sample_details_01[, c(1,3,4)] %>% 
  full_join(knn_imputed_data, by = "sample") %>% 
  select(- "sample")

## 2.2. for metabo data ----------------------------------------

load(metabo_data)

Data_metabo <- sample_details_01[, c(1,3,4)] %>% 
  full_join(metabo_data, by = "sample") %>% 
  select(- "sample")

# 3. Function for calculating the variance removal effect -----------------------------

residualFactor <- function(Data, sample_data, Factor, save_residual_model_feeding_regime){
  # initiating linear model
  response_colname <- colnames(sample_data)[-1]
  treatment_colname <-  c("Feeding_Regime", "Slaughter_Condition")
  treatment_colname_ref <- c("Control", "NoStress")
  
  # creating reference/control
  for (i in seq_along(treatment_colname)){
    Data[[treatment_colname[i]]] <- factor(Data[[treatment_colname[i]]])
    Data[[treatment_colname[i]]] <- relevel(Data[[treatment_colname[i]]], ref = treatment_colname_ref[i])
  }
  
  # producing the residual (data without variation of the feeding_regime)
  
  residual_model_feeding_regime <- data.frame()
  
  for (i in 1: length(response_colname)){
    # Correct formula construction
    formula <- as.formula(glue("{response_colname[i]} ~ {Factor}"))
    model_feeding_regime <- lm(formula, data = Data)
    
    residual_model_feeding_regime[1:40,response_colname[i]] <-  resid(model_feeding_regime)
    
  }
  
  residual_model_feeding_regime <- residual_model_feeding_regime %>% 
    mutate(sample = sample_data$sample) %>% 
    select(sample, everything())
  
  write.csv(residual_model_feeding_regime, file = save_residual_model_feeding_regime, row.names = FALSE)
}


# 3.  Removing variation due to factor "Feeding_Regime" from knn_imputed data -------------------
Factor <- "Feeding_Regime"
residualFactor(Data_knn_imputed, knn_imputed_data, Factor, save_residual_model_feeding_regime_knn_imputed)


# 4. Removing variation due to factor "Feeding_Regime" from metabo data -------------------
Factor <- "Feeding_Regime"
residualFactor(Data_metabo, metabo_data, Factor, save_residual_model_feeding_regime_metabo)


