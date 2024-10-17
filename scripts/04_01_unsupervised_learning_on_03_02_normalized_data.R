# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: ""data/12_anova_preFold_result_table_post_06.RData"
# Tutorials:

# 1. Libraries ------------------------------------------------------------

# Basic libraries
pacman::p_load(tidyverse, 
               ggplot2, 
               plotly, # for interactive plotting
               glue
               )


#for PCA
# basic stat package and ggplot is enough, loading plot need to be created

#t-SNE
pacman::p_load(Rtsne, tsne)

#UMAP
pacman::p_load(umap) 

# Basic libraries
pacman::p_load(tidylog)


source("scripts/00_path_variables.R")
source("scripts/functions/getDimensionReduction.R")

# 2. Importing data --------------------------------------------

## 2.1. Importing metabo_normalized data ----------------------------------
load(metabo_data)
df <-  metabo_data %>% 
  column_to_rownames("sample")

## 2.2. Importing sample labels ------------------------------
load(sample_details_01)
labels <- sample_details_01 %>% 
  column_to_rownames("sample")


# 3. Dimensionality Reduction ---------------------------------------------

## 3.1. PCA -----------------

### 3.1.1. Defining PCs to plot ---------------------

x = 1 # PC1
y = 2 # PC2
model = "pca"

### 3.1.2. Model building -------
pca_model <- prcomp(df, scale. = TRUE) # my data is already pareto scaled, but still this scaling is to get unit variance

### 3.1.3. Extracting PCs -----
pca_data <- as.data.frame(pca_model$x)

### 3.1.4. Calculating explained variance ------
explained_variance <- pca_model$sdev^2 / sum(pca_model$sdev^2)
explained_variance <- round(explained_variance * 100, 2)

### 3.1.5. Preparing Scree Plot -----

scree_plot_data <- data.frame(PC_component = paste0('PC ', 1: nrow(metabo_data)), explained_variance = cumsum(explained_variance))
components_to_show = 10

scree_plot <- generate_pca_scree_plot(scree_plot_data, components_to_show)
scree_plot

### 3.1.6. Preparing Score Plots -------

for (label in colnames(labels)) {
  assign(paste0("fig_pca_scorePlot_", label), generate_dimension_reduction_plot("pca",x, y, pca_data, label, labels, explained_variance = explained_variance))
}

fig_pca_scorePlot_label
fig_pca_scorePlot_Slaughter_Condition
fig_pca_scorePlot_Feeding_Regime





## 3.2. t-SNE (https://kb.10xgenomics.com/hc/en-us/articles/217265066-What-is-t-Distributed-Stochastic-Neighbor-Embedding-t-SNE) ----------

### 3.2.1. Defining parameters to plot -----
x = 1
y = 2 
perplexity = 20
model = "tsne"

### 3.2.2. Model building -------  
set.seed(0)
tsne_model <- tsne(df, initial_dims = 2)
# tsne_model <- Rtsne(df, perplexity = perplexity, pca_scale = TRUE)

### 3.2.3. Extracting tSNE components -----
tsne_data <- data.frame(tsne_model)
# tsne_data <- as.data.frame(tsne_model$Y)

rownames(tsne_data) <- rownames(df)


### 3.2.4. preparing tSNE plots -----

for (label in colnames(labels)) {
  assign(paste0("fig_tsne_Plot_", label), generate_dimension_reduction_plot(model, x, y, tsne_data, label, labels, explained_variance = NULL))
}

fig_tsne_Plot_label
fig_tsne_Plot_Slaughter_Condition
fig_tsne_Plot_Feeding_Regime


## 3.3. UMAP -----------
### 3.3.1. Defining parameters to plot -----
x = 1
y = 2 
model = "umap"

### 3.2.2. Model building -------  
set.seed(0)
umap_model <- umap(df, n_components = 2, random_state = 15)

### 3.2.3. Extracting UMAP components -----
umap_data <- data.frame(umap_model[["layout"]]) 

### 3.2.4. preparing UMAP plots -----

for (label in colnames(labels)) {
  assign(paste0("fig_umap_Plot_", label), generate_dimension_reduction_plot(model, x, y, umap_data, label, labels, explained_variance = NULL))
}

fig_umap_Plot_label
fig_umap_Plot_Slaughter_Condition
fig_umap_Plot_Feeding_Regime



# 4. Saving images as ggplot ----------------------------------------------

items <- ls()

for (item in items){
  # print(item)
  if(grepl("fig", item)){
    print(glue("plots/{item}.png"))
    ggsave(get(item), file = glue("plots/04_01_{item}.png"), width = 8, height = 8, units = "in", dpi = 600)
  }
}


