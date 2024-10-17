# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: 
# Input
# csv_data_for_metabo <- "data/03_01_03_data_for_metabo_analysis_normalization_after_03_01.csv"
# output
# save_metabo_data <- "data/03_02_01_metabo_data_for_statistical_analysis_and_enrichment_post_03_02.RData" 

# The data is being normalized also sample details were populated with all the treatment names, eg. Horn, Stress

# 1. Libraries ------------------------------------------------------------
# BiocManager::install(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea"))


# # impute: impute.knn(); pcaMethods: different types of PCA implementation, https://www.bioconductor.org/packages/release/bioc/vignettes/pcaMethods/inst/doc/pcaMethods.pdf;
# # globaltest: factorial ANOVA (maybe); GlobalAncova: Global test for groups of variables via model comparisons; multtest: Resampling-based multiple hypothesis testing; fgsea: goana and kegga for GO ad KEGG analysis
# 
# pacman::p_load(devtools)
# devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

pacman::p_load(MetaboAnalystR, tidyverse, tidylog)
source("scripts/00_path_variables.R")


# 2. Pre-processing with MetaboAnalyst ------------------------------------

mSet<-
  # Concentration and Statistics Package of MetaboAnalyst
  InitDataObjects("conc", "stat", FALSE) %>%  
  # rowu =  row-unpaired
  Read.TextData(.,csv_data_for_metabo, "rowu", "disc") %>%
  # Checking data 
  SanityCheckData(.) %>% 
  ### They are converting the data some way I am not able to catch, for kNN Imputation
  #"A total of 40208 (59%) missing values were detected." 
  #Remove features with more than 10% missing values
  RemoveMissingPercent(., percent=0.1) %>%
  # Imputing missing variables with knn from "impute" package
  ImputeMissingVar(., method="knn_var")%>% 
  # Post imputation they are performing the Fold Change calculation
  # SanityCheckData(.) %>% #  FilterVariable(mSet, "F", 10, "iqr", 0, "mean", 0) %>%
  PreparePrenormData(.) %>%
  Normalization(., "NULL", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)

metabo_data = mSet[["dataSet"]][["norm"]] %>% 
  rownames_to_column(var = "sample")

save(metabo_data, file = save_metabo_data)


# 3. Density plot with limma to check the distribution -----------------
grDevices::windows(16, 10)
par(mfrow = c(1,2))
gene_wise_density_plot <- limma::plotDensities(metabo_data[, -1], 
                     legend = FALSE,
                     main = "gene wise density plot")
gene_wise_density_plot

sample_wise_density_plot <- limma::plotDensities(t(metabo_data[, -1]), 
                     legend = FALSE, 
                     main = "sample wise density plot") 
sample_wise_density_plot




####################################################################################### sample details should be automatically there

# 
# # 3. Cleaning brackets of sample_details$sample and converting to more descriptive labels --------
# load(sample_details)
# 
# ## 3.1. Cleaning brackets of sample_details$sample --------------
# source("scripts/functions/removeBracket.R") # from "sample_details$sample" has bracket in it, 
# # where as Metaboanalyst remove it by default. So, converting E-PP-BIOM-8712(31) -> E-PP-BIOM-871231
# 
# sample_details$sample <- removeBracket(sample_details$sample)
# 
# ## 3.2. Converting to more descriptive labels --------------
# source("scripts/functions/cleanLabels.R") # this to convert "NoStress_NCM -> NoStress"
# 
# # Define the label mappings for each class
# label_mappings <- list(
#   "Horn" = list(
#     list(pattern = "NC", new_label = "NoHorn"),
#     list(pattern = ".*", new_label = "Horn")  # Default label if no other pattern matches
#   ),
#   "Stress" = list(
#     list(pattern = "NoStress", new_label = "NoStress"),
#     list(pattern = ".*", new_label = "Stress")
#   ),
#   "Rearing" = list(
#     list(pattern = "NM", new_label = "NoMixed"),
#     list(pattern = ".*", new_label = "Mixed")
#   )
# )
# 
# # Adding the clean labels to the "sample_details" data.frame
# sample_details_01 = add_column(sample_details, 
#                                RC =  as.vector(cleanLabels(sample_details$label, "Rearing", label_mappings)), 
#                                SC =  as.vector(cleanLabels(sample_details$label, "Stress", label_mappings)),
#                                Horns =  as.vector(cleanLabels(sample_details$label, "Horn", label_mappings)),.after = 2)
# 
# 
# # 4. Merging sample_details with the metabo_data --------------------------
# 
# metabo_data <- sample_details %>% 
#   full_join(metabo_data, by= join_by("sample"))
# 
# # 3. Saving "metabo_data" for further analysis ----------------------------
# 
# # save(metabo_data, file = "data/10_metabo_data_for_statistical_analysis_and_enrichment_post_03_02.RData")
# # save(sample_details_01, file = "data/11_sample_details_01_with_all_the_trt_names_post_03_02.RData")
