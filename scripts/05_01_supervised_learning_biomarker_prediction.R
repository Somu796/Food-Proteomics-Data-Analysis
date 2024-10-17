#### I have remove variation caused by other variables than my slaughter condition....https://www.perplexity.ai/search/i-have-two-factors-can-i-remov-55_eWFQLTWab_iE2Zh.2BA#0 
# Ensemble Feature Selection with EFS in R

# examples
# https://cran.r-project.org/web/packages/EFS/EFS.pdf
# https://d-nb.info/1140232363/34
# https://biodatamining.biomedcentral.com/articles/10.1186/s13040-016-0114-4#citeas


# Sys.setenv(BIOCONDUCTOR_CONFIG_FILE="http://bioconductor.org/config.yaml")

# 1. install and import the package ---------------------------------------
# Genze N, Neumann U (2017). EFS: Tool for Ensemble Feature Selection. R package version 1.0.3, https://CRAN.R-project.org/package=EFS

# install.packages("EFS")
pacman::p_load(tidyverse, glue, janitor , EFS, tidylog)
source("scripts/00_path_variables.R")


# 2. Importing Data ----------------------------------------------

# A. General Approach
## sample details
load(sample_details_01)

## sample data
sample_data <- qs::qread("data_proc.qs")%>% 
  rownames_to_column(var = "sample")

## Data merging
Data <- sample_details_01[, c(1,3,4)] %>% 
  full_join(sample_data, by = "sample") %>% 
  select(- "sample")

# ## Preparing data for classification (in general situation)
# Data <- Data %>% 
#   select(- Feeding_Regime)


# B. Manual Approach for  lipivimus data
## Data from Variation removal due to "Feeding_Regime"
residual_model_feeding_regime <- read.csv(save_residual_model_feeding_regime)

## Preparing data for classification
data_efs <- residual_model_feeding_regime %>% 
  mutate(Slaughter_Condition = Data$Slaughter_Condition) %>% 
  select("Slaughter_Condition", everything()) %>% 
  select(- "sample")




data_efs$Slaughter_Condition <- ifelse(Data$Slaughter_Condition == "Stress", 1, 0)


# 4. Building the Ensemble Feature Selection ------------------------------

# Start feature selection
efs <- ensemble_fs(data = data_efs, classnumber = 1,
                   NA_threshold = 0.2, cor_threshold = 0.7,
                   runs = 100) #, selection = rep(TRUE, 8)
view(as.data.frame(efs))

# saving the result in an excel sheet
efs_result <- c()
efs_result <- efs %>% 
  t() %>% 
  as.data.frame() 

efs_result <- efs_result %>% 
  mutate(total = rowSums(.)) %>% 
  arrange(desc(total)) %>% 
  rownames_to_column(var = "gene_name")

# Save the file
# efs_result_dest <- "data/05_01_Ensemble_Feature_selection_output.csv"
write.csv(efs_result, file = efs_result_dest,row.names = FALSE)


# Create a cumulative barplot based on the output from EFS
barplot_fs("docs/slaughter_condition", efs)

# # Create a ROC Curve based on the output from efs
# eval_tests <- efs_eval(data = Data, efs_table = efs,
#                        file_name = "docs/Horn",
#                        classnumber = 1, NA_threshold = 0.2,
#                        logreg = TRUE,
#                        permutation = TRUE, p_num = 100,
#                        variances = TRUE, jaccard = TRUE,
#                        bs_num = 100, bs_percentage = 0.9)


