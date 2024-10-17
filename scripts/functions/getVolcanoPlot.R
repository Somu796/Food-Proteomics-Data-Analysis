library(ggplot2)
library(extrafont)
loadfonts(device = "win")
library(ggrepel)
library(glue)
pacman::p_load(ggpp)

# VolcanoPlot(Factor, FOldChangeData, mrlimit = 1.2, pValue_Limit= 0.01, pAdjust = TRUE)

getVolcanoPlot = function(Factor, meanFactorLevel.1, meanFactorLevel.2, FOldChangeData, mrlimit = 2, pValue_Limit= 0.05, pAdjust = FALSE, pAdjustMethod = "fdr", legend_title = NULL){
  #on-off proteins(only present in one group
  #Factor eg. Horns, RC, SC
  
  differential_data <- data.frame(gene_name = FOldChangeData$gene_name , 
                                  meanRatio = FOldChangeData$MeanRatio, 
                                  pValue = FOldChangeData[[Factor]], 
                                  differential_group = "No") # forming a new differential data for plotting
  # Adding a row to inform Up regulated or Down regulated
  differential_data$differential_group[(differential_data$meanRatio> mrlimit) & differential_data$pValue< pValue_Limit] = "Up"
  differential_data$differential_group[(differential_data$meanRatio< 1/mrlimit) & differential_data$pValue< pValue_Limit] = "Down"
  
  # Making a factor of the differential group to assign color to it
  differential_data = differential_data %>%
    mutate(differential_group = factor(differential_group, levels = c("Up", "No", "Down")))
  
  # Assigning the color
  regulation_color <- setNames(c("#FF7F0E", "grey", "#1F77B4"), # #bb0c00
                               levels(differential_data$differential_group))
  
  # getting the plots for plotting
  differential_data$delabel <- ifelse(differential_data$differential_group == "Up" | differential_data$differential_group == "Down", 
                                      differential_data$gene_name, NA)
  # Total up, down and not sig entries 
  total_Up = sum(differential_data$differential_group == "Up")
  print(total_Up)
  total_Down = sum(differential_data$differential_group == "Down")
  print(total_Down)
  total_Not_Sig = sum(is.na(differential_data$delabel))
  print(total_Not_Sig)
  
  # p-Adj Values Calculation
  if (pAdjust == TRUE){
    differential_data$fdr.p.value <- p.adjust(differential_data$pValue, method = pAdjustMethod)
  }
  
  # Plotting the result
  
  
  volcano_plot = ggplot(data = differential_data, 
                        mapping = aes(log2(meanRatio), 
                                      -log10(pValue), 
                                      col= differential_group #, 
                                      # label = delabel
                                      )) +
    geom_point(alpha = 5/10, 
               size = 2) +
    geom_vline(xintercept = c(-log2(mrlimit), 
                              log2(mrlimit)), 
               col = "black", 
               linetype = 'dashed') +
    geom_hline(yintercept = -log10(pValue_Limit), 
               col = "black", 
               linetype = 'dashed') + 
    # to overcome the text overlap
    # geom_text_repel(aes(family = "serif"),
                    # max.overlaps = Inf
                    # position =
                    #   position_nudge_to(x = 2.3),
                    # min.segment.length = 0,
                    # segment.color = "black",
                    # arrow = arrow(length = unit(0.015, "npc")),
                    # direction = "y"
                    # ) +
    labs(color = legend_title, 
         x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = meanFactorLevel.1, C = meanFactorLevel.2)) ,
         y = expression("-log"[10]*"(Pvalue)"),
         caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit} \n Up = {total_Up}; Down = {total_Down}; Not sig. = {total_Not_Sig}")) +
    scale_color_manual(values = regulation_color, # to set the colours of our variable  
                       labels = c("Up", "Not sig.", "Down")) +
    theme_minimal() + 
    theme(text = element_text(family = "serif"),
          axis.line = element_line(),
          axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
          axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
          plot.title = element_text(hjust = 0.5)
    )
  
  return(volcano_plot)
  
}


#FF7F0E



















# VolcanoPlot = function(Factor, meanRatio, pValue, mrlimit = 2, pValue_Limit= 0.05, pAdjust = FALSE, pAdjustMethod = "fdr"){
#   #on-off proteins(only present in one group
#   #Factor eg. Horns, RC, SC
#   
#   differential.mr <- (meanRatio > mrlimit | meanRatio < 1/mrlimit)
#   differential.p <- (pValue < pValue_Limit)
#   differential <- (differential.mr & differential.p)
#   which.differential <- which(differential)
#   # Differential Proteins
#   
#   if (length(which.differential) == 0) print("There are NO Differential Proteins") else cat(sprintf("There are %d Differential Proteins", length(which.differential)))
#   # p-Value Adjustment
#   if (pAdjust == TRUE){
#     fdr.p.value <- p.adjust(pValue, method = pAdjustMethod)
#   }
#   
#   #Plotting 1
#   plot.data <- cbind(log.mr = log2(meanRatio), log.p = -log10(pValue))
#   
#   ggplot(data = plot.data, mapping = aes(log.mr, log.p)) +
#     geom_point(alpha = 5/10) +
#     geom_point()
#   
#   
#   plot.data <- cbind(log.mr = log2(meanRatio), log.p = -log10(pValue))
#   main <- paste("Volcano plot for", Factor)
#   xlab <- "log2(fold change)"
#   ylab <- "- log10 (p-value)"
#   plot(plot.data, col = "grey", main = main, xlab = xlab, ylab = ylab, pch = 19)
#   abline(h = -log10(0.05), lty = 2, col = "grey30",
#          v = log2(c(mrlimit, 1/mrlimit)))
#   points(plot.data[differential,], col = "black", pch = 19)
#   #text(perfor, policy, canton, pos=1)
#   
#   # post multiple testing correction
#   best.proteins = which(fdr.p.value < 0.05)
#   
#   if (length(best.proteins) == 0) print("There are NO best proteins") else cat(sprintf("There are %d Best Proteins", length(best.proteins)))
#   points(plot.data[best.proteins,], col = "darkorange2", pch = 19)
# }
# 




