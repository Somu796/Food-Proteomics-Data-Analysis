# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: 
# Input
# csv_data_for_metabo <- "data/03_01_03_data_for_metabo_analysis_normalization_after_03_01.csv"
# output
# save_metabo_data <- "data/03_02_01_metabo_data_for_statistical_analysis_and_enrichment_post_03_02.RData" 

# The data is being normalized also sample details were populated with all the treatment names, eg. Horn, Stress
# the mean and median I was calclulating in a wrong way as it was not n=mean and median under each population i.e., type of sample

# 1. Libraries ------------------------------------------------------------
pacman::p_load(tidyverse, limma, rstatix, car, glue, moments, tidylog)
source("scripts/00_path_variables.R")
load(metabo_data)
load(sample_details_01)

# colnames(metabo_data) <- gsub("-", "_",colnames(metabo_data)) # 'SERPINA3-2' was giving error

# Data for using in formulas and groupby
Data <- sample_details_01 %>% 
  full_join(metabo_data, by = "sample")

# Data for using in complete numeric operation like shapiro wilk test
Data_Numeric <-  metabo_data %>% 
  column_to_rownames("sample")

# 2. Checking overall Normality with Density Plot ---------------------------------------

## 2.1. Checking all together ----------------------------------------------

grDevices::windows(16, 10)
par(mfrow = c(1,2))
gene_wise_density_plot <- limma::plotDensities(Data_Numeric, 
                                               legend = FALSE,
                                               main = "gene wise density plot")
gene_wise_density_plot

sample_wise_density_plot <- limma::plotDensities(t(Data_Numeric), 
                                                 legend = FALSE, 
                                                 main = "sample wise density plot") 
sample_wise_density_plot


## 2.2. Checking one by one
sample_no <- 30
gene_no <- 50

grDevices::windows(16, 10)
par(mfrow = c(1,2))
gene_wise_density_plot <- limma::plotDensities(Data_Numeric[, gene_no], 
                                               legend = FALSE,
                                               main = glue("gene wise density plot (gene name : {colnames(Data_Numeric)[gene_no]})"))
gene_wise_density_plot

sample_wise_density_plot <- limma::plotDensities(t(Data_Numeric)[sample_no, ], 
                                                 legend = FALSE, 
                                                 main = glue("sample wise density plot (sample number : {rownames(Data_Numeric)[sample_no]})")) 
sample_wise_density_plot

# par(mfrow = c(3,5))
# hist(metabo_data[,1])
# hist(metabo_data[,2])
# hist(metabo_data[,3])
# hist(metabo_data[,4])
# hist(metabo_data[,5])
# hist(metabo_data[,6])
# hist(metabo_data[,7])


# 3. Mean == Median check (Without grouping sample from different population doesn't make sense, the data is not coming from same population) -----------------------------------------------------------------------

summary_data <-Data %>% 
  group_by(label) %>% 
  get_summary_stats(show = c("mean", "median"))


summary_data$diff_mean_median <- summary_data$mean - summary_data$median

colnames(summary_data)[2] <- "Accession"

# Distribution of difference between mean and median to see if the difference is a lot for assumption checking ANOVA
groups <- unique(summary_data$label)
grDevices::windows(16, 10)
par(mfrow = c(2, 4))

for (i in 1: length(groups)){
  x <- summary_data %>% 
    filter(summary_data$label == groups[i]) %>% 
    pull(diff_mean_median)
  
  hist(x, main = paste("Histogram of" , groups[i]), xlab = "diff value (mean - median)")
  # title(groups[i])
  }

# 4.Calculating skewness (using moments package of R)  ---------------------------------------------------------------------
skewness_result <- Data %>% 
  group_by(label) %>% 
  reframe(across(
    .cols = colnames(Data)[5: length(colnames(Data))], 
    .fns = skewness,
    .names = "{col}")) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = colnames(Data)[5] : colnames(Data)[length(colnames(Data))],
               names_to = "Accession",
               values_to = "skewness") 
  

# 5. Shapiro Wilk Test ----------------------------------------------------

shapiro_result <- Data %>%
  group_by(label) %>%
  reframe(across(
    .cols = colnames(Data)[5:length(colnames(Data))], # Apply from 5th column onwards
    .fns = ~ shapiro.test(.x)$p.value,                # Extract only p-value
    .names = "{col}"                             # Naming format for p-value columns
  )) %>%
  as.data.frame()%>% 
  pivot_longer(cols = colnames(Data)[5] : colnames(Data)[length(colnames(Data))],
               names_to = "Accession",
               values_to = "shapiro_p_value_before") 

# shapiro_data <- Data %>%
#   group_by(label) %>%
#   reframe(across(
#     .cols = colnames(Data)[5:length(colnames(Data))],
#     .fns = ~ list(shapiro.test(.x)), # Apply Shapiro-Wilk test and store results in list
#     .names = "{col}"
#   )) %>%
#   as.data.frame()



# 6. Homogenity of Variance (Levene Test) -----------------------------------------------
# this is an over all test for all the groups

levene_result <- c()
for (i in colnames(Data_Numeric)){
  levene_result[[i]] = leveneTest(as.formula(glue("{i} ~ Feeding_Regime*Slaughter_Condition")), data = Data)$`Pr(>F)`[1]
  
}

leveneTest(as.formula(glue("{i} ~ Feeding_Regime*Slaughter_Condition")), data = Data)


levene_result <- levene_result %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Accession")

colnames(levene_result)[2] <- "levene_test"


# 7. Merging all the data -------------------------------

assumption_checking <- summary_data %>% 
  full_join(skewness_result, by = c("label", "Accession")) %>% 
  full_join(shapiro_result, by = c("label", "Accession")) %>% 
  # full_join(shapiro_test_after, by = "Accession") %>%
  full_join(levene_result, by = "Accession") 




# 8. Be awake about this variables ----------------------------

assumption_checking %>% 
  filter(shapiro_p_value_before <0.001) %>% 
  view()

# 9. Saving Assumption Checking Data --------------------------
save(assumption_checking, file = "data/12_01_assumption_cheking_pre_anova.RData")