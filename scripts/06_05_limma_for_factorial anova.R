## Documentation
# https://rformassspectrometry.github.io/QFeatures/articles/Processing.html#imputation-1 :: QFetaure
# https://rformassspectrometry.github.io/book/sec-quant.html :: R for Mass Spectrometry
# https://uclouvain-cbio.github.io/WSBIM2122/sec-prot.html :: Omics Data Analysis
# https://lgatto.github.io/MSnbase/articles/v01-MSnbase-demo.html :: MSnbase, Normalization

# Conclusion: came really bad for my data only 17 differential expression, for mixed it is coming a lot, verify with cpact data following  
# https://uclouvain-cbio.github.io/WSBIM2122/sec-prot.html

# 1. Loading Libraries ----------------------------------------------------

# loading libraries
pacman::p_load(tidyverse, glue, SummarizedExperiment, QFeatures, limma, msdata, factoextra, naniar)
source("scripts/00_path_variables.R")

# 2. Loading Data ---------------------------------------------------------

## 2.1. row Details and assay ----------------
load(enriched_gene_data_01)


rownames(enriched_gene_data_01) <- enriched_gene_data_01$Accession

enriched_gene_data_01 <- enriched_gene_data_01 %>% 
  mutate(across(c(colnames(enriched_gene_data_01)[8:47]), ~as.numeric(.)))

## 2.2. column details --------------
load(sample_details_01)

sample_details_01 <-  sample_details_01 %>%
  dplyr::rename(quantCols = sample)  # quantCols is considered to be key between assay and column data


sample_details_01$Feeding_Regime <- factor(sample_details_01$Feeding_Regime) %>% 
  relevel(ref = "Control")

sample_details_01$Slaughter_Condition <- factor(sample_details_01$Slaughter_Condition) %>% 
  relevel(ref = "NoStress")

# Combinations exists
table(sample_details_01$Feeding_Regime, sample_details_01$Slaughter_Condition)

# 3. Reading My Data with Qfeatures ---------------------------------------

## 3.1. Using Shiny App -----------------

# QFeaturesGUI::importQFeatures()

## 3.2. Using Code -------------------

data <- c()
data <- readQFeatures(enriched_gene_data_01, 
                      quantCols = 8:47, 
                      fnames = "Accession", 
                      colData = sample_details_01, 
                      name = "proteins") #, fnames = "Accession"

data[[1]]

# adding colData in summarized Experiment because QFeature only add to itself
colData(data[["proteins"]]) <- colData(data) # ths line adds up colData in colData names field

# Data Information 
data[["proteins"]]

# row data information
rowData(data[[1]])

# col data information
colData(data[["proteins"]]) #summarizedExperiment
colData(data) # QFeature bject

# assay data
assay(data[["proteins"]]) %>% 
  as.data.frame() %>% 
  view()

# 4. Missing Values -------------------------------------------------------
## 4.1. Visualizing Misising Values --------------------
data[["proteins"]] <- zeroIsNA(data[["proteins"]]) # take summarized Experiment Object
NA_data <- nNA(data[["proteins"]])  

# heatmap for visualizing the data  
# heatmap with overall data
vis_miss(as.data.frame(t(assay(data[["proteins"]]))), sort_miss = TRUE, show_perc_col = FALSE) +
  labs(x = "Genes") +
  theme(axis.text.x = element_blank())

# heatmap with only have missing data
col_with_missingval <- as.data.frame(t(assay(data[["proteins"]]))) %>%
  select(where(~ any(is.na(.))))

vis_miss(col_with_missingval,   sort_miss = TRUE, show_perc = FALSE) +
  labs(x = "Genes") +
  theme(axis.text.x = element_text(face="bold", size = 7))

# barplot
barplot(NA_data$nNAcols$pNA)

# table of the datahttp://127.0.0.1:9055/graphics/plot_zoom_png?width=1920&height=1017
table(NA_data$nNArows$nNA)

# # sorting gene_name vs number of sample has missing values  
# NA_data$nNArows %>% 
#   as.data.frame() %>% 
#   view()
# 
# total_miss_col <- NA_data$nNArows %>% 
#   as.data.frame() %>% 
#   filter(pNA>=0.1) %>% 
#   nrow()
# print(glue("Total missing column {total_miss_col}"))
# 
# 
# total_col <- NA_data$nNArows %>%
#   nrow()
# print(glue("Total column {total_col}"))
# 
# left_col <- total_col - total_miss_col
# print(glue("Total left column {left_col}"))

## 4.2. imputing missing values -----------------
# remove rows that have more than 10% missing data
data <- addAssay(data,
                 filterNA(data[["proteins"]], pNA = 0.1),
                 name = "NAfiltered_proteins")


NA_data <- nNA(data[["NAfiltered_proteins"]]) 
NA_data
# data[["proteins"]] <- filterNA(data[["proteins"]], pNA = 0.1) # Alternative

nNA(data[["proteins"]])


# 5. Log Transform ----------------------------------------------
# this is just to verify the plots
new_data <- addAssay(data,
                     logTransform(data[["NAfiltered_proteins"]]),
                     name = "log_NAfiltered_proteins")

# 6. Imputing Missing Values ----------------------------------------------------------------------

## 6.1. Plotting how imputed data follow the distribution with the original data ------------------------- 
# plotting with normal data
par(mfrow = c(1, 2))
cls <- c("black", "red", "blue", "steelblue", "orange")
plot(density(na.omit(assay(data[["NAfiltered_proteins"]]))), col = cls[1])

lines(density(assay(impute(data[["NAfiltered_proteins"]], method = "knn"))), col = cls[2])
lines(density((assay(impute(data[["NAfiltered_proteins"]], method = "zero")))), col = cls[3])
lines(density(assay(impute(data[["NAfiltered_proteins"]], method = "MinDet"))), col = cls[4])
lines(density(assay(impute(data[["NAfiltered_proteins"]], method = "bpca"))), col = cls[5])
legend("topright", legend = c("orig", "knn", "zero", "MinDet", "bpca"),
       col = cls, lwd = 2, bty = "n")

# plotting with log transform data for better visualization
cls <- c("black", "red", "blue", "steelblue", "orange")
plot(density(na.omit(assay(new_data[["log_NAfiltered_proteins"]]))), col = cls[1])

lines(density(assay(impute(new_data[["log_NAfiltered_proteins"]], method = "knn"))), col = cls[2])
lines(density((assay(impute(new_data[["log_NAfiltered_proteins"]], method = "zero")))), col = cls[3])
lines(density(assay(impute(new_data[["log_NAfiltered_proteins"]], method = "MinDet"))), col = cls[4])
lines(density(assay(impute(new_data[["log_NAfiltered_proteins"]], method = "bpca"))), col = cls[5]) # This won't work if some column has more than one 
legend("topright", legend = c("orig", "knn", "zero", "MinDet", "bpca"),
       col = cls, lwd = 2, bty = "n")

## 6.2. Imputing data with kNN algorithm -------------------------------
# imputation with kNN
data <- addAssay(data,
                 impute(data[["NAfiltered_proteins"]], method = "knn"),
                 name = "knnimp_NAfiltered_proteins"
)

# plot(data)

# An instance of class QFeatures containing 3 assays:
# [1] ST_proteomics: SummarizedExperiment with 750 rows and 79 columns 
# [2] log_ST_proteomics: SummarizedExperiment with 750 rows and 79 columns 
# [3] imputedAssay: SummarizedExperiment with 750 rows and 79 columns

# verifying imputation is done
nNA(data[["knnimp_NAfiltered_proteins"]])$nNA$nNA

# 7. Log transformation after imputation  ---------------------------------------------------

data <- addAssay(data,
                 logTransform(data[["knnimp_NAfiltered_proteins"]]),
                 name = "log_knnimp_NAfiltered_proteins")

# 8. Normalisation --------------------------------------------------------

data <- addAssay(data,
                 normalize(data[["log_knnimp_NAfiltered_proteins"]],
                           method = "center.median"),
                 name = "norm_log_knnimp_NAfiltered_proteins")

# for sample
grDevices::windows(16, 10)
par(mfrow = c(2, 2))
limma::plotDensities(assay(data[["proteins"]]), legend = FALSE)
limma::plotDensities(assay(data[["knnimp_NAfiltered_proteins"]]), legend = FALSE)
limma::plotDensities(assay(data[["log_knnimp_NAfiltered_proteins"]]), legend = FALSE)
limma::plotDensities(assay(data[["norm_log_knnimp_NAfiltered_proteins"]]), legend = FALSE)

# 9. PCA -------------

# label
pca_prot_label <-
  data[["norm_log_knnimp_NAfiltered_proteins"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp() %>%
  fviz_pca_ind(habillage = data$label,
               title = "PCA")
pca_prot_label

# All
pca_prot_all <-
  data[["norm_log_knnimp_NAfiltered_proteins"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp() 

pca_prot_all$x%>%
  as.data.frame() %>%
  ggplot(aes(x = .[, 1], y = .[, 2], color = data$Slaughter_Condition, shape = data$Feeding_Regime)) +
  geom_point() + theme_minimal()

# Slaughter COndition
pca_prot_slaughter <-
  data[["norm_log_knnimp_NAfiltered_proteins"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp() %>%
  fviz_pca_ind(habillage = data$Slaughter_Condition,
               title = "Proteins (median aggregation)")

pca_prot_slaughter

# 10. Statistical Analysis ------------------------------------------------

## 10.1. Extracting data  from the QFeature ---------------------------
prots <- getWithColData(data, "norm_log_knnimp_NAfiltered_proteins")

## 10.2. Performing Limma ---------------------------

### 10.2.1. Classic Interaction model with sum to zero parametrization  ---------------------------

#### 10.2.1.1. creating design matrix ---------------------------
contrasts(prots$Slaughter_Condition) <- contr.sum(2)
contrasts(prots$Feeding_Regime) <- contr.sum(4)

design <- model.matrix(~prots$Slaughter_Condition * prots$Feeding_Regime)

# removed because that is coming as NA 
design <- design[,-8]

# fitting with limma
fit <- lmFit(assay(prots), design)
fit$coefficients # print to see which is producing NA
fit_Bayes <- eBayes(fit)

top_proteins <- topTable(fit_Bayes, n=Inf, adjust = "BH")
#
#### 10.2.1.2. Getting the  main effects and contrasts ----------------------------------

top_proteins_stress <- topTable(fit_Bayes,  n=Inf, coef = "prots$Slaughter_Condition1", adjust = "none")
# volcanoplot(top_proteins_stress)

# plot(top_proteins_stress$logFC, -log10(top_proteins_stress$P.Value))
# 
# volcanoplot(fit_Bayes, coef = "prots$Slaughter_Condition1", style = "p-value", highlight = 0,  hl.col="blue",
#             xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

top_proteins_diet <- topTable(fit_Bayes,  n=Inf, coef = c("prots$Feeding_Regime1", "prots$Feeding_Regime2", "prots$Feeding_Regime3"), adjust = "none")


### 10.2.1. Analysis as a single factor  ---------------------------

# making a single factor
group <- paste(prots$Slaughter_Condition, prots$Feeding_Regime, sep = "_")

# preparing the design
design <- model.matrix(~0+group)

# modifying the design colnames
colnames(design) <- gsub("group","", colnames(design))

# 
contrasts_group <- makeContrasts(
  Stress_vs_NoStress <- (Stress_Lipid +Stress_Lipid_VitE +Stress_Lipid_VitE_PlantExt)/3 - (NoStress_Control + NoStress_Lipid + NoStress_Lipid_VitE + NoStress_Lipid_VitE_PlantExt)/4,
  
  Lipid_vs_Control_NoStress <- NoStress_Lipid - NoStress_Control,
  Lipid_VitE_vs_Control_NoStress <- NoStress_Lipid_VitE - NoStress_Control,
  Lipid_VitE_PlantExt_vs_Control_NoStress <- NoStress_Lipid_VitE_PlantExt - NoStress_Control,
  Lipid_VitE_vs_Lipid_NoStress <- NoStress_Lipid_VitE - NoStress_Lipid,
  Lipid_VitE_PlantExt_vs_Lipid_NoStress <- NoStress_Lipid_VitE_PlantExt - NoStress_Lipid,
  Lipid_VitE_PlantExt_vs_Lipid_VitE_NoStress <- NoStress_Lipid_VitE_PlantExt - NoStress_Lipid_VitE,
  
  Lipid_VitE_vs_Lipid_Stress <- Stress_Lipid_VitE - Stress_Lipid,
  Lipid_VitE_PlantExt_vs_Lipid_Stress <- Stress_Lipid_VitE_PlantExt - Stress_Lipid,
  Lipid_VitE_PlantExt_vs_Lipid_VitE_Stress <- Stress_Lipid_VitE_PlantExt - Stress_Lipid_VitE,
  
  
  levels = colnames(design)
)

model <- lmFit(assay(prots), design)

model_contr <-  contrasts.fit(model, contrasts_group)
model_contr_bayes <- eBayes(model_contr)

top2_proteins <- topTable(model_contr_bayes, n=Inf, adjust = "none")

top2_proteins_stress <- topTable(model_contr_bayes, coef = 1, n=Inf, adjust = "none")


### In case of unbalaced data the values doesn't match.
top_proteins_stress %>% 
  filter(P.Value <0.05) %>% 
  nrow()

top2_proteins_stress %>% 
  filter(P.Value <0.05) %>% 
  nrow()
