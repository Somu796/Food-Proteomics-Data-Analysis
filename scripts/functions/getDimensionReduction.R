
# 0. general Dimension reduction Plot (for PCA, t-SNE, UMAP a general scatter plot function) ---------------------

generate_dimension_reduction_plot <- function(model, x, y, model_output_data, label, labels, explained_variance = NULL) {
  
  # pca data and corresponding "sample" name
  model_output_data <- model_output_data  %>%
    rownames_to_column("sample")
  
  # group label and corresponding "sample" name
  # label <- "Feeding_Regime"
  group_names <- labels %>% 
    select(all_of(label))%>% 
    rownames_to_column("sample")
  
  colnames(group_names)[2] <- "group" 
  
  # joining pca data and group label by "sample"
  model_output_data <- model_output_data %>% 
    full_join(group_names, by = "sample") %>% 
    column_to_rownames("sample")
  
  
  print(model_output_data)
  
    # add_column(group = Data[[label]])
  
  p <-
    ggplot(model_output_data,
           aes(x = model_output_data[, x], y = model_output_data[, y], color = group)) +
    geom_point(size = 4, alpha = 5 / 10) +  # Marker size
    geom_hline(yintercept = 0, linetype="dashed", alpha = 3/10) +
    geom_vline(xintercept = 0, linetype="dashed", alpha = 3/10) +
    geom_polygon(stat = "ellipse",
                 aes(fill = group),
                 alpha = 0.001) +
    theme_bw(base_family = "serif")
  
  if (tolower(model) == "pca") {
    p <- p +
      labs(
        x = glue("Principal Component {x} ({explained_variance[x]}%)"),
        y = glue("Principal Component {y} ({explained_variance[y]}%)")
      )
    
  } else if (tolower(model) == "tsne"){
    p <- p  +
      labs(
        x = glue("tSNE {x}"),
        y = glue("tSNE {y}")
      )
  } else if (tolower(model) == "umap"){
    p <- p  +
      labs(
        x = glue("UMAP {x}"),
        y = glue("UMAP {y}")
      )
  } else{
    p <- p  +
      labs(
        x = glue("Component {x}"),
        y = glue("Component {y}")
      )
  }
  
  # show(p)
  # fig_scorePlot <- ggplotly(p)
  # fig_scorePlot
  
  return(p)
}

# 1. PCA ------------------------------------------------------------------

## 1.1. PCA Scree Plot -------------------

generate_pca_scree_plot <- function(scree_plot_data, n){
  scree_plot <- ggplot(scree_plot_data[1:n,], mapping = aes(x = reorder(PC_component, explained_variance), y = explained_variance)) +
    geom_bar(stat = "identity", fill="#56B4E9", colour="black") +
    geom_line(aes(x = 1:n, y = explained_variance), linewidth = 1) +
    geom_text(aes(label = round(explained_variance, 2)), vjust = -0.5, size = 3)+
    labs(x = "Principal Componenets",
         y = "Explained Variance (%)") +
    theme_bw(base_family = "serif")
  return(scree_plot)
}


## 1.2. PCA Score Plot -------------------

generate_pca_score_plot <- function(x, y, pca_data, label, explained_variance) {
  pca_data <- pca_data %>% 
    add_column(group = Data[[label]])
  
  p <- ggplot(pca_data, aes(x = pca_data[,x], y = pca_data[,y], color = group)) +
    geom_point(size = 4, alpha = 5/10) +  # Marker size
    geom_polygon(stat = "ellipse", 
                 aes(fill = group), 
                 alpha = 0.001)+ 
    labs(x = glue("Principal Component {x} ({explained_variance[x]}%)"),
         y = glue("Principal Component {y} ({explained_variance[y]}%)")) +
    theme_bw()
  
  # show(p)
  # fig_scorePlot <- ggplotly(p)
  # fig_scorePlot
  
  return(p)
}

# Testing
# generate_pca_score_plot(x, y, pca_data, label, explained_variance)

# 2. t-SNE --------------------

generate_tSNE_plot <- function(x, y, tSNE_data, label) {
  
  tsne_data <- tsne_data %>% 
    add_column(group = Data[[label]])
  
  p <- ggplot(tsne_data, aes(x = tsne_data[,x], y = tsne_data[,y], color = group)) +
    geom_point(size = 4, alpha = 5/10) +  # Marker size
    geom_polygon(stat = "ellipse", 
                 aes(fill = group), 
                 alpha = 0.001)+ 
    labs(x = glue("tSNE {x}"),
         y = glue("tSNE {y}")) +
    theme_bw()
  
  # show(p)
  # fig_scorePlot <- ggplotly(p)
  # fig_scorePlot
  
  return(p)
}



