# Description
# To RUN: This is an independent script. This has to be run manually. 
# saving_data: yes
# data_name:  "data/06_02_02_factorial_anova_posthoc_result.csv"
# This calculates the factorial anova and group by mean for fold change calculation
# limma tutroials: https://www.youtube.com/live/HWOllPTpnq0?feature=shared
# creating design matrices: https://f1000research.com/articles/9-1444
# original paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402510/
# linear models concept paper: https://doi.org/10.2202/1544-6115.1027

# WARNING: genenames might not b unique better to use the acccession insted of Accession

# To automate: 

# 1. Libraries ------------------------------------------------------------
pacman::p_load(qs, tidyverse, janitor, broom, limma, car, emmeans, glue, tidylog) # rstatix,
source("scripts/00_path_variables.R")


# 2. importing Data -------------------------------------------------------
load(metabo_data)
load(sample_details_01)
load(gene_data_02)

# colnames(metabo_data) <- gsub("-", "_",colnames(metabo_data)) # 'MT-...' was giving error

# Data for using in formulas and groupby
Data <- sample_details_01 %>% 
  full_join(metabo_data, by = "sample")

# write.csv(Data, file= "data/minitab/normalized_data.csv")

# 3. Input for the Statistical Test ---------------------------------------
response_colname <- colnames(metabo_data)[-1]
treatment_colname <-  c("Feeding_Regime", "Slaughter_Condition")
treatment_colname_ref <- c("Lipid", "NoStress")

post_hoc_columns <- "Feeding_Regime" # because it has 5 output
p_limit <- 1 # Consider all the p-Vaues

# 4. Defining reference factors and creating sum contrsts (https://stackoverflow.com/a/68742165/25712971) -----------------

## 4.1. creating reference/control

for (i in seq_along(treatment_colname)){
  Data[[treatment_colname[i]]] <- factor(Data[[treatment_colname[i]]])
  Data[[treatment_colname[i]]] <- relevel(Data[[treatment_colname[i]]], ref = treatment_colname_ref[i])
}

# 5. Performing Anova and Post-hoc tests ---------------------------------------------------

## 5.1. Performing car::Anova ---------------------------------------------
factorial_anova_result <- data.frame()
  # Accession = rep(NA, ncol(metabo_data)-1), 
  #                                     p_Feeding_Regime = rep(NA, ncol(metabo_data)-1), 
  #                                     p_Slaughter_Condition = rep(NA, ncol(metabo_data)-1))


# To get the main Effects
for (i in 1: length(response_colname)){
  # Correct formula construction
  formula <- as.formula(glue("{response_colname[i]} ~ {paste(treatment_colname, collapse = '*')}"))
  model <- lm(formula, data = Data)
  # summary(model)
  anova.result_III <- tidy(car::Anova(model, type = 3))
  if (anova.result_III[anova.result_III$term == "Feeding_Regime:Slaughter_Condition",]$p.value >0.05){
    anova.result <- tidy(car::Anova(model, type = 2))
  }else{
    anova.result <- anova.result_III
  }
  
  # tidy(anova.result)
  print(response_colname[i])
  factorial_anova_result[i, "Accession"] <- response_colname[i]
  
  for (term_id in 1: length(anova.result$term)){
    factorial_anova_result[[anova.result$term[term_id]]][i] <- anova.result$p.value[term_id]  
  }
}

factorial_anova_result <- gene_data_02 %>% 
  select(Accession, gene_name, Protein_name) %>% 
  right_join(factorial_anova_result, by = "Accession")


write.csv(factorial_anova_result, file= save_factorial_anova_result, row.names = FALSE)

# anova.result <- tidy(car::Anova(model, type = 3))

## 5.2. Performing Post-Hoc on Feeding_Regime --------------------------------
# To get the post-hoc result

proteins_with_sig_feeding_regime <- factorial_anova_result[factorial_anova_result$Feeding_Regime < 0.05, ]$Accession
factorial_anova_posthoc_result <- data.frame()

comparison_formula <-  as.formula("~ Feeding_Regime|Slaughter_Condition")

for (i in 1: length(proteins_with_sig_feeding_regime)){
  formula <- as.formula(glue("{response_colname[i]} ~ {paste(treatment_colname, collapse = '*')}"))
  model <- lm(formula, data = Data)
  posthoc_model <- emmeans(model, comparison_formula)
  posthoc_result <-  summary(pairs(posthoc_model))
  posthoc_result[["Treatment"]] <- paste(posthoc_result[,2], posthoc_result[,1], sep= " : ")
  
  print(response_colname[i])
  factorial_anova_posthoc_result[i, "Accession"] <- response_colname[i]
  
  for (treatment_id in 1: length(posthoc_result[["Treatment"]])){
    # print(posthoc_result$p.value[treatment_id])
    factorial_anova_posthoc_result[[posthoc_result$Treatment[treatment_id]]][i] <- posthoc_result$p.value[treatment_id]
  }
}

factorial_anova_posthoc_result <- gene_data_02 %>% 
  select(Accession, gene_name, Protein_name) %>% 
  right_join(factorial_anova_posthoc_result, by = "Accession")

write.csv(factorial_anova_posthoc_result, file= save_factorial_anova_posthoc_result, row.names = FALSE)
