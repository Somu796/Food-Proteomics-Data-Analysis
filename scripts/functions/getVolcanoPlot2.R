# This is Volcano Plot with EFS

library(ggplot2)
library(extrafont)
loadfonts(device = "win")
library(ggrepel)
library(glue)
pacman::p_load(ggpp)



getVolcanoPlotEFS <- function(Factor, meanFactorLevel.1, meanFactorLevel.2, FOldChangeData, mrlimit = 2, pValue_Limit= 0.05, pAdjust = FALSE, pAdjustMethod = "fdr", legend_title = NULL, additional_data = NULL){
  if(is.null(additional_data)){
    print("additional_data is not provided")
    break
  }else{
    differential_data <- data.frame(gene_name = FOldChangeData$gene_name , 
                                    meanRatio = FOldChangeData$MeanRatio,
                                    pValue = FOldChangeData[[Factor]],
                                    differential_group = "No") # forming a new differential data for plotting
    
    # Color
    differential_data$EFS <- "No"
    differential_data[(differential_data$meanRatio> mrlimit) & differential_data$pValue< pValue_Limit, ]$EFS = "Sig"
    differential_data[(differential_data$meanRatio< 1/mrlimit) & differential_data$pValue< pValue_Limit, ]$EFS = "Sig"
    differential_data[differential_data$gene_name %in% additional_data$gene_name, ]$EFS <- "Yes"
     
    differential_data$EFS <- factor(differential_data$EFS, levels = c("Yes", "Sig", "No"))
    regulation_color <- setNames(c("red", "black", "grey"),
                                levels(differential_data$differential_group))
    # Text 
    differential_data$delabel <- ifelse(differential_data$EFS == "Yes",
                                        differential_data$gene_name, NA)
    # Plotting
    volcano_plot = ggplot(data = differential_data,
                          mapping = aes(log2(meanRatio),
                                        -log10(pValue),
                                        label = delabel)) +
      # point
      geom_point(
        aes(colour = EFS), 
        # shape= 21,
        size = 2,
        alpha = 5/10
        )  +
      # line
      geom_vline(xintercept = c(-log2(mrlimit),
                                log2(mrlimit)),
                 col = "black",
                 linetype = 'dashed') +
      geom_hline(yintercept = -log10(pValue_Limit),
                 col = "black",
                 linetype = 'dashed')+
      #scale_fill_manual(values = rep("white", 3))   +
      # text: to overcome the text overlap
      geom_text_repel(max.overlaps = Inf, 
                      color =  "red",
                      position = 
                        position_nudge_to(x = 2.3),
                      min.segment.length = 0,
                      segment.color = "black",
                      arrow = arrow(length = unit(0.015, "npc")),
                      direction = "y") +
      labs(#color = legend_title,
        x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = meanFactorLevel.1, C = meanFactorLevel.2)) ,
        y = expression("-log"[10]*"(Pvalue)"),
        caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit}")) +
      scale_color_manual(values = regulation_color) + 
      
      
      theme_minimal() +
      theme(text = element_text(family = "Times New Roman"),
            axis.line = element_line(),
            axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
            axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
            plot.title = element_text(hjust = 0.5)
      )

    return(volcano_plot)
    
  }
  
  }

# getVolcanoPlotEFS(Factor, meanFactorLevel.1, meanFactorLevel.2 ,FOldChangeData, mrlimit, pValue_Limit, pAdjust = TRUE, legend_title = "Horn Status", additional_data = filtered_data)
# 
# regulation_color <- setNames(c("red", "black", "grey"),
#                              levels(differential_data$differential_group))
# 
# 
# 
# 
# 
# volcano_plot = ggplot(data = differential_data,
#                       mapping = aes(log2(meanRatio),
#                                     -log10(pValue),
#                                     label = delabel)) +
#   # point
#   geom_point(aes(colour = EFS), shape=21, size = 2)  +
#   # line
#   geom_vline(xintercept = c(-log2(mrlimit),
#                             log2(mrlimit)),
#              col = "black",
#              linetype = 'dashed') +
#   geom_hline(yintercept = -log10(pValue_Limit),
#              col = "black",
#              linetype = 'dashed')+
#   #scale_fill_manual(values = rep("white", 3))   +
#   # text: to overcome the text overlap
#   geom_text_repel(max.overlaps = Inf, 
#                   color =  "red",
#                   point.padding = 0.2, 
#                   nudge_x = .15,
#                   nudge_y = .5,
#                   segment.curvature = -1e-20,
#                   arrow = arrow(length = unit(0.015, "npc"))) +
#   labs(#color = legend_title,
#        x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = meanFactorLevel.1, C = meanFactorLevel.2)) ,
#        y = expression("-log"[10]*"(Pvalue)"),
#        caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit}")) +
#   scale_color_manual(values = regulation_color) + 
#   
#   
#   theme_minimal() +
#   theme(text = element_text(family = "Times New Roman"),
#         axis.line = element_line(),
#         axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
#         axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
#         plot.title = element_text(hjust = 0.5)
#   )
# 
# 
# #### Try more
# 
# ggplot(data = differential_data,
#        mapping = aes(log2(meanRatio),
#                      -log10(pValue),
#                      label = delabel)) +
#   # point
#   geom_point(aes(colour = EFS), shape=21, size = 2)  +
#   # line
#   geom_vline(xintercept = c(-log2(mrlimit),
#                             log2(mrlimit)),
#              col = "black",
#              linetype = 'dashed') +
#   geom_hline(yintercept = -log10(pValue_Limit),
#              col = "black",
#              linetype = 'dashed')+
#   #scale_fill_manual(values = rep("white", 3))   +
#   # text: to overcome the text overlap
#   geom_text_repel(max.overlaps = Inf, 
#                   color =  "red",
#                   position = 
#                     position_nudge_to(x = 2.3),
#                   min.segment.length = 0,
#                   segment.color = "black",
#                   arrow = arrow(length = unit(0.015, "npc")),
#                   direction = "y") +
#   labs(#color = legend_title,
#     x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = meanFactorLevel.1, C = meanFactorLevel.2)) ,
#     y = expression("-log"[10]*"(Pvalue)"),
#     caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit}")) +
#   scale_color_manual(values = regulation_color) + 
#   
#   
#   theme_minimal() +
#   theme(text = element_text(family = "Times New Roman"),
#         axis.line = element_line(),
#         axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
#         axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
#         plot.title = element_text(hjust = 0.5)
#   )
# 
# ## Try Again
# 
# ggplot(data = differential_data,
#        mapping = aes(log2(meanRatio),
#                      -log10(pValue),
#                      label = delabel)) +
#   # point
#   geom_point(aes(colour = EFS), shape=21, size = 2)  +
#   # line
#   geom_vline(xintercept = c(-log2(mrlimit),
#                             log2(mrlimit)),
#              col = "black",
#              linetype = 'dashed') +
#   geom_hline(yintercept = -log10(pValue_Limit),
#              col = "black",
#              linetype = 'dashed')+
#   #scale_fill_manual(values = rep("white", 3))   +
#   # text: to overcome the text overlap
#   geom_text_repel( 
#     color = "red",
#     verbose = TRUE,
#     seed = 123,
#     max.time = 1,
#     max.iter = Inf,
#     size = 3) +
#   labs(#color = legend_title,
#     x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = meanFactorLevel.1, C = meanFactorLevel.2)) ,
#     y = expression("-log"[10]*"(Pvalue)"),
#     caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit}")) +
#   scale_color_manual(values = regulation_color) + 
#   
#   
#   theme_minimal() +
#   theme(text = element_text(family = "Times New Roman"),
#         axis.line = element_line(),
#         axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
#         axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
#         plot.title = element_text(hjust = 0.5)
#   )
# 
# 
# 
# filtered_data$gene_name
# differential_data <- data.frame(gene_name = FOldChangeData$gene_name , 
#                                 meanRatio = FOldChangeData$MeanRatio,
#                                 pValue = FOldChangeData[[Factor]],
#                                 differential_group = "No") # forming a new differential data for plotting
# 
# 
# 
# volcano_plot = ggplot(data = differential_data,
#                       mapping = aes(log2(meanRatio),
#                                     -log10(pValue),
#                                     colour = EFS,
#                                     fill = "black", 
#                                     label = delabel)) +
#   geom_point()
# 
# 
# df = data.frame(x = 1:10, y = 11:20, col = factor(rep(c("a", "b"), 5)))
# 
# ggplot(df, aes(x, y))+ geom_point(aes(colour=factor(col)), shape=21, size = 2) + 
#   scale_fill_manual(values=c("white", "white")) + 
#   scale_colour_manual(values=c("red", "black")) +
# 
# 
# 
# 
# 
# differential_data$meanRatio
# # Color
# differential_data$EFS <- "No"
# differential_data[(differential_data$meanRatio> mrlimit) & differential_data$pValue< pValue_Limit, ]$EFS = "Sig"
# differential_data[(differential_data$meanRatio< 1/mrlimit) & differential_data$pValue< pValue_Limit, ]$EFS = "Sig"
# differential_data[differential_data$gene_name %in% filtered_data$gene_name]$EFS <- "Yes"
# 
# differential_data$EFS <- factor(differential_data$EFS, levels = c("Yes", "Sig", "No"))
# regulation_color <- setNames(c("red", "black", "grey"),
#                              levels(differential_data$differential_group))
# # Text 
# differential_data$delabel <- ifelse(differential_data$EFS == "Yes" | differential_data$EFS == "Sig",
#                                     differential_data$gene_name, NA)
# # Plotting
# volcano_plot = ggplot(data = differential_data,
#                       mapping = aes(log2(meanRatio),
#                                     -log10(pValue),
#                                     colour = EFS,
#                                     label = delabel))
# 
# 
# 
# # library(ggplot2)
# # library(extrafont)
# # loadfonts(device = "win")
# # library(ggrepel)
# # library(glue)
# # 
# # # VolcanoPlot(Factor, FOldChangeData, mrlimit = 1.2, pValue_Limit= 0.01, pAdjust = TRUE)
# # 
# # getVolcanoPlot = function(Factor, meanFactorLevel.1, meanFactorLevel.2, FOldChangeData, mrlimit = 2, pValue_Limit= 0.05, pAdjust = FALSE, pAdjustMethod = "fdr", legend_title = NULL, additional_data = NULL){
# #   #on-off proteins(only present in one group
# #   #Factor eg. Horns, RC, SC
# #   
# #   differential_data <- data.frame(gene_name = FOldChangeData$gene_name , 
# #                                   meanRatio = FOldChangeData$MeanRatio, 
# #                                   pValue = FOldChangeData[[Factor]], 
# #                                   differential_group = "No") # forming a new differential data for plotting
# #   # Adding a row to inform Up regulated or Down regulated
# #   differential_data$differential_group[(differential_data$meanRatio> mrlimit) & differential_data$pValue< pValue_Limit] = "Up"
# #   differential_data$differential_group[(differential_data$meanRatio< 1/mrlimit) & differential_data$pValue< pValue_Limit] = "Down"
# #   
# #   # Making a factor of the differential group to assign color to it
# #   differential_data = differential_data %>%
# #     mutate(differential_group = factor(differential_group, levels = c("Up", "No", "Down")))
# #   
# #   # Assigning the color
# #   regulation_color <- setNames(c("#FF7F0E", "grey", "#1F77B4"), # #bb0c00
# #                                levels(differential_data$differential_group))
# #   
# #   # getting the plots for plotting
# #   differential_data$delabel <- ifelse(differential_data$differential_group == "Up" | differential_data$differential_group == "Down", 
# #                                       differential_data$gene_name, NA)
# #   # Total up, down and not sig entries 
# #   total_Up = sum(differential_data$differential_group == "Up")
# #   print(total_Up)
# #   total_Down = sum(differential_data$differential_group == "Down")
# #   print(total_Down)
# #   total_Not_Sig = sum(is.na(differential_data$delabel))
# #   print(total_Not_Sig)
# # 
# #   
# #   # plotting addditional data 
# #   
# #   
# #   # Plotting the result
# #   volcano_plot = ggplot(data = differential_data, 
# #                         mapping = aes(log2(meanRatio), 
# #                                       -log10(pValue),
# #                                       colour = differential_group)) +
# #     geom_point(alpha = 5/10, 
# #                size = 2)
# #   if(!is.null(additional_data)){
# #     volcano_plot = volcano_plot +
# #       geom_point(data = additional_data, 
# #                  aes(log2(MeanRatio),
# #                      -log10(additional_data[[Factor]]),
# #                      fill = "yellow",
# #                      label = delabel), size = 3)
# #   }
# #   
# #   volcano_plot <- volcano_plot +
# #     geom_vline(xintercept = c(-log2(mrlimit), 
# #                               log2(mrlimit)), 
# #                col = "black", 
# #                linetype = 'dashed') +
# #     geom_hline(yintercept = -log10(pValue_Limit), 
# #                col = "black", 
# #                linetype = 'dashed') + 
# #     # to overcome the text overlap
# #     geom_text_repel(max.overlaps = Inf) +  
# #     labs(color = legend_title, 
# #          x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = meanFactorLevel.1, C = meanFactorLevel.2)) ,
# #          y = expression("-log"[10]*"(Pvalue)"),
# #          caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit} \n Up = {total_Up}; Down = {total_Down}; Not sig. = {total_Not_Sig}")) +
# #     scale_color_manual(values = regulation_color, # to set the colours of our variable  
# #                        labels = c("Up", "Not sig.", "Down")) +
# #     theme_minimal() + 
# #     theme(text = element_text(family = "Times New Roman"),
# #           axis.line = element_line(),
# #           axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
# #           axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
# #           plot.title = element_text(hjust = 0.5)
# #     )
# # 
# #   return(volcano_plot)
# #   
# # }
# # 
# # 
# # #FF7F0E
# # 
# # 
# # important_genes_parts = genes <- c(
# #   "ALDH7A1", "ETAA1", "DDX1", "COPS6", "CCT3", "ALDH6A1", "IDH2", "FGG",
# #   "USP9X", "HSP90B1", "PHB2", "ACTC1", "GPT", "COX1", "DDX39B",
# #   "BTF3", "SAMD15", "SURF4", "EHD2", "SCRN3"
# # )
# # 
# # 
# # filtered_data <- FOldChangeData[FOldChangeData$gene_name %in% important_genes_parts,]
# # 
# # filtered_data$delabel <- filtered_data$gene_name
# # 
# # getVolcanoPlot(Factor, meanFactorLevel.1, meanFactorLevel.2 ,FOldChangeData, mrlimit, pValue_Limit, pAdjust = TRUE, legend_title = "Horn Status", additional_data = filtered_data)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # VolcanoPlot = function(Factor, meanRatio, pValue, mrlimit = 2, pValue_Limit= 0.05, pAdjust = FALSE, pAdjustMethod = "fdr"){
# # #   #on-off proteins(only present in one group
# # #   #Factor eg. Horns, RC, SC
# # #   
# # #   differential.mr <- (meanRatio > mrlimit | meanRatio < 1/mrlimit)
# # #   differential.p <- (pValue < pValue_Limit)
# # #   differential <- (differential.mr & differential.p)
# # #   which.differential <- which(differential)
# # #   # Differential Proteins
# # #   
# # #   if (length(which.differential) == 0) print("There are NO Differential Proteins") else cat(sprintf("There are %d Differential Proteins", length(which.differential)))
# # #   # p-Value Adjustment
# # #   if (pAdjust == TRUE){
# # #     fdr.p.value <- p.adjust(pValue, method = pAdjustMethod)
# # #   }
# # #   
# # #   #Plotting 1
# # #   plot.data <- cbind(log.mr = log2(meanRatio), log.p = -log10(pValue))
# # #   
# # #   ggplot(data = plot.data, mapping = aes(log.mr, log.p)) +
# # #     geom_point(alpha = 5/10) +
# # #     geom_point()
# # #   
# # #   
# # #   plot.data <- cbind(log.mr = log2(meanRatio), log.p = -log10(pValue))
# # #   main <- paste("Volcano plot for", Factor)
# # #   xlab <- "log2(fold change)"
# # #   ylab <- "- log10 (p-value)"
# # #   plot(plot.data, col = "grey", main = main, xlab = xlab, ylab = ylab, pch = 19)
# # #   abline(h = -log10(0.05), lty = 2, col = "grey30",
# # #          v = log2(c(mrlimit, 1/mrlimit)))
# # #   points(plot.data[differential,], col = "black", pch = 19)
# # #   #text(perfor, policy, canton, pos=1)
# # #   
# # #   # post multiple testing correction
# # #   best.proteins = which(fdr.p.value < 0.05)
# # #   
# # #   if (length(best.proteins) == 0) print("There are NO best proteins") else cat(sprintf("There are %d Best Proteins", length(best.proteins)))
# # #   points(plot.data[best.proteins,], col = "darkorange2", pch = 19)
# # # }
# # # 
# # 
# # 
# # 
# # 
