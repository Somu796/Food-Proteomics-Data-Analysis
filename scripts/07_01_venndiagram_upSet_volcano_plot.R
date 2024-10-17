# Description
# To RUN: This is an independent script.
# saving_data: YES
# data_name: 
# This is to calculate the Upset Plot and volcanoplot
# Reading Materials:
# venn diagram: https://www.datanovia.com/en/blog/beautiful-ggplot-venn-diagram-with-r/
# upSet Plot: https://krassowski.github.io/complex-upset/articles/Examples_R.html


# 1. Loading Libraries ----------------------------------------------------
pacman::p_load(xlsx, tidyverse, ggplot2, ggvenn, UpSetR, ggrepel, pheatmap, glue, grDevices, svglite,tidylog)
# library(extrafont)
# loadfonts(device = "win")

# 2. Loading the data -----------------------------------------------------

source("scripts/00_path_variables.R")
factorial_anova_result <- read.csv(save_factorial_anova_result) # for venn diagram
factorial_anova_posthoc_result <- read.csv(save_factorial_anova_posthoc_result) # for upset plot
load(fold_change_result)
load(gene_data_02)
load(sample_details_01)
load(metabo_data)

# 3. Venn Diagram ---------------------------------------------------------

# venn_diagram_data <- list(
#   `Feeding Regime` = factorial_anova_result[factorial_anova_result$Feeding_Regime < 0.05, "protein"],
#   `Slaughter Condition` = factorial_anova_result[factorial_anova_result$Slaughter_Condition < 0.05, "protein"],
#   `Feeding Regime * Slaughter Condition` = factorial_anova_result[factorial_anova_result$Feeding_Regime.Slaughter_Condition < 0.05, "protein"]
# )

source(changingVennFont)

venn_diagram_data <- tibble(
  protein = factorial_anova_result$protein,
  `Feeding Regime` = factorial_anova_result$Feeding_Regime < 0.05,
  `Slaughter Condition` = factorial_anova_result$Slaughter_Condition < 0.05,
  `Feeding Regime : Slaughter Condition` = factorial_anova_result$Feeding_Regime.Slaughter_Condition < 0.05
)

## 3.1. venn diagram with 2 main effect and one simple main effect
venn_plot_3 <- ggplot(venn_diagram_data, aes(A = `Feeding Regime`,
                                   B = `Slaughter Condition`,
                                   C = `Feeding Regime : Slaughter Condition`)) +
  geom_venn(text_size = 5, set_name_size = 7, show_outside = "none", fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))  +
  theme_void() +
  theme(text = element_text(family = "Comic Sans MS"))+
  coord_fixed()

png(save_venn_diagram_interaction, units="in", width=8, height=8, res=600)
par(family="serif")
# venn_plot_3
changingVennFont(venn_plot_3, font = "serif")
dev.off()

## 3.2. venn diagram with only 2 main effects
venn_plot_2 <- ggplot(venn_diagram_data, aes(A = `Feeding Regime`,
                                   B = `Slaughter Condition`)) +
  geom_venn(text_size = 5, set_name_size = 7, show_outside = "none",  fill_color = c("#0073C2FF", "#EFC000FF"),
            stroke_size = 1)  +
  theme_void() +
  theme(text = element_text(family = "Comic Sans MS"))+
  coord_fixed()


png(save_venn_diagram_wo_interaction, units="in", width=8, height=8, res=600)
par(family="serif")
changingVennFont(venn_plot_2, font = "serif")
dev.off()


# 4. Upset Plot with Post-hoc Results --------------------------------------

pValue_Limit= 0.05

## 4.1. Separating the No Stress and Stress df separately -------------------------------
nostress_colnames <- grepl("NoStress", colnames(factorial_anova_posthoc_result))
stress_colnames <- setdiff(colnames(factorial_anova_posthoc_result), colnames(factorial_anova_posthoc_result)[nostress_colnames])
nostress_colnames[1] <- TRUE

## 4.2. Modifying the column names -----------------------

### 4.2.1. replacing "..." with "-"  ----------------------
nostress_df <- factorial_anova_posthoc_result[,nostress_colnames]
colnames(nostress_df) <- str_replace(colnames(nostress_df),"^NoStress\\.\\.\\.", "") %>% 
  str_replace(., "\\.\\.\\."," \\- ")

stress_df <- factorial_anova_posthoc_result[,stress_colnames]
colnames(stress_df) <- str_replace(colnames(stress_df),"^Stress\\.\\.\\.", "") %>% 
  str_replace(., "\\.\\.\\."," \\- ")

### 4.2.2.  swapping "stressA-stressB" to "stressB-stressA" ---------------------------

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

## 4.3. Making a column for upset plot ----------------------------

stress_df[is.na(stress_df)] <- 2 # replacing NA with 2

# getting the column to plot upset
for (i in 1: nrow(stress_df)){
  stress_df$posthoc[i] <- str_replace_all(
    paste(
      unlist(
        colnames(stress_df)[stress_df[i,] <0.05]
      ), collapse = "&"
    ), "_", "+")
}


nostress_df[is.na(nostress_df)] <- 2 # replacing NA with 2 

# getting the column to plot upset
for (i in 1: nrow(nostress_df)){
  nostress_df$posthoc[i] <- str_replace_all(
    paste(
      unlist(
        colnames(nostress_df)[nostress_df[i,] <0.05]
        ), collapse = "&")
    , "_", "+")
}


## 4.4. Upset plot with R ------------------------------------------------

### 4.4.1. Removing set size from upSetR plot--------

# skip_set_size_plot <- function(ups) {
#   main <- ups$Main_bar
#   ## find panel grob
#   panel_idx <- grep("panel", main$layout$name, fixed = TRUE)
#   ## find text grob
#   text_idx <- which(
#     vapply(main$grobs[[panel_idx]]$children, 
#            \(x) inherits(x, "text"), 
#            logical(1)))
#   tG <- main$grobs[[panel_idx]]$children[[text_idx]]
#   tG$label <- paste0(tG$label, " (",
#                      scales::label_percent(0.1)(as.numeric(tG$label) /
#                                                   sum(as.numeric(tG$label))),
#                      ")")
#   main$grobs[[panel_idx]]$children[[text_idx]] <- tG
#   grid.newpage()
#   grid.draw(arrangeGrob(main, ups$Matrix, heights = ups$mb.ratio))
# }
# 

### 4.4.1. Upset plot for the Stress -------------------------------------------
upSet_data_stress <- stress_df %>% 
  filter(posthoc != "") %>% 
  select(posthoc) %>% 
  table()

#Plot
# windows(16,10)
upset_plot_diet_givenstress <- upset(fromExpression(upSet_data_stress),
                         # nsets = 3,
                         mb.ratio = c(0.6, 0.4),
                         number.angles = 0,
                         # text.scale = 1.1,
                         point.size = 2.8,
                         line.size = 1,
                         # sets.bar.color = c("HDC"="#9ad8a1",
                         #                   "CE"="#f5a3b7",
                         #                   "Raw" = "#9ce9f8",
                         #                   "HPP" = "#9cadea",
                         #                   "UAE" = "#f794f2")
                         # order.by = "degree",
                         # decreasing = FALSE,
                         # set_size.numbers_size = TRUE,
                         text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2.5) # last: numbers above bar, 2nd last set names
)
# upset_plot_diet_givenstress
# skip_set_size_plot(upset_plot_diet_givenstress)

png(save_upset_plot_diet_givenstress, units="in", width=14, height=8, res=600)
par(family="serif")
upset_plot_diet_givenstress
dev.off()

### 4.4.2. Upset plot for NoStress -----------------
upSet_data_nostress <- nostress_df %>% 
  filter(posthoc != "") %>% 
  select(posthoc) %>% 
  table()


# windows(16,10)
upset_plot_diet_givennostress <- upset(fromExpression(upSet_data_nostress),
                                     # nsets = 3,
                                     mb.ratio = c(0.6, 0.4),
                                     number.angles = 0,
                                     # text.scale = 1.1,
                                     point.size = 2.8,
                                     line.size = 1,
                                     # sets.bar.color = c("HDC"="#9ad8a1",
                                     #                   "CE"="#f5a3b7",
                                     #                   "Raw" = "#9ce9f8",
                                     #                   "HPP" = "#9cadea",
                                     #                   "UAE" = "#f794f2")
                                     # order.by = "degree",
                                     # decreasing = FALSE,
                                     # set_size.numbers_size = TRUE,
                                     # sets = NA,
                                     # sets.bar.color = NA,
                                     # sets.x.label = "",
                                     # scale.sets = "",
                                     text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2.5) # last: numbers above bar, 2nd last set names
)
# upset_plot_diet_givennostress


png(save_upset_plot_diet_givennostress, units="in", width=14, height=8, res=1200)
par(family="serif")
upset_plot_diet_givennostress
dev.off()




# 5. Heatmap for significant proteins -------------------------------------
### Moved it to 07_03 page

# 6. Volcano Plot ---------------

source(getFoldChange)
source(getVolcanoPlot)


factorial_anova_result <- gene_data_02[, 1:2] %>% 
  right_join(factorial_anova_result, by = "Accession") %>% 
  full_join(fold_change_result, by = "Accession")


mrlimit = 1.2
pValue_Limit= 0.05

Factor = "Slaughter_Condition"
meanFactorLevel.1 = "Stress"
meanFactorLevel.2 = "NoStress"

columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
FOldChangeData = c()
FOldChangeData = getFoldChange(factorial_anova_result, columns)

# write.csv(FOldChangeData, file = glue("data/07_01_{Factor}_anova_foldchange_result_post_07.csv"), row.names = FALSE)

sc_volcano = getVolcanoPlot(Factor, meanFactorLevel.1, meanFactorLevel.2, FOldChangeData, mrlimit, pValue_Limit, pAdjust = TRUE, legend_title = "Slaughter Condition")
# windows(10,10)
# sc_volcano

ggsave(sc_volcano, file = save_volcanoplot_slaughter_condition, dpi = 600)

# 7. upStream and DownStream Regulator ----------------------
# Up and Down Regulated Genes
upregulated_genes <- FOldChangeData %>% 
  filter(Slaughter_Condition<0.05, log2FC>log2(1.2))

downregulated_gene <- FOldChangeData %>% 
  filter(Slaughter_Condition<0.05,log2FC< -log2(1.2))

write.xlsx(upregulated_genes, file = save_up_and_down_regulated_genes, sheetName = "upregulated_genes", append = TRUE)
write.xlsx(downregulated_gene, file = save_up_and_down_regulated_genes, sheetName = "downregulated_gene", append = TRUE)

