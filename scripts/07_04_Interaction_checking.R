pacman::p_load(tidyverse)
source("scripts/00_path_variables.R")

df <- read.csv("data/06_02_01_factorial_anova_result.csv")

colnames <- df %>% 
  filter(Slaughter_Condition<0.05 & Feeding_Regime.Slaughter_Condition<0.05) %>% 
  pull(Accession) #%>% 
  # paste(collapse = ", ")

# colnames(df)

load(metabo_data)
load(sample_details_01)
# load(gene_data_02)

# colnames(metabo_data) <- gsub("-", "_",colnames(metabo_data)) # 'MT-...' was giving error

# Data for using in formulas and groupby
Data <- sample_details_01 %>% 
  full_join(metabo_data, by = "sample")

md <-  Data %>% 
  select("Slaughter_Condition", "Feeding_Regime", all_of(colnames))

# md$new <-  paste(md$Slaughter_Condition, md$Feeding_Regime, )

# write.csv(md, "data/minitab/02_normalized_data.csv", row.names = FALSE)


  
md$new <- paste(md$Slaughter_Condition, md$Feeding_Regime, sep = ".")
response_colname <-  colnames(md)[c(-1,-2, -ncol(md))]

# 
# i <- 1
# m <- lm(as.formula(glue("{response_colname[i]} ~ Feeding_Regime*Slaughter_Condition")), data = md)
# posthoc_model <- emmeans(m, 
#                          as.formula("pairwise~ Feeding_Regime|Slaughter_Condition"))
# posthoc_result <-  summary(pairs(posthoc_model))
# posthoc_result
# 
# emmip(m, Slaughter_Condition~Feeding_Regime)
# plot(posthoc_result, xlab = response_colname[i])
# 
# 
# 
# gb <-  md %>% 
#   group_by()


### Write in Console
# installed.packages(c("readxl", "multcompView","dplyr", "ggplot2"))

### Run the code below
# Calling needed libraries
# library(readxl)
library(ggplot2)
library(multcompView)
# library(dplyr)
library(car)

i <-  1

md$Slaughter_Condition <- factor(md$Slaughter_Condition, levels = c("NoStress", "Stress"))
md$Feeding_Regime <- factor(md$Feeding_Regime, levels = unique(md$Feeding_Regime))

unique(md$Feeding_Regime)

for (i in 1: length(response_colname)){

  # # ANOVA One-Way
  anova.result = aov(as.formula(glue("{response_colname[i]} ~ Slaughter_Condition* Feeding_Regime")), data = md)
    # car::Anova(lm(as.formula(glue("{response_colname[i]} ~ new")), data = md), type = "III")
  # summary(anova.result)
  
  # Summarizing the Data in mean and Standard deviation
  Data_summary =  md %>% 
    group_by(Slaughter_Condition, Feeding_Regime) %>%
    reframe(across(
      .cols = i, # Select columns 3 to the last column dynamically
      .fns = list(
        mean = mean,
        sd = sd),
      .names = "{.fn}" # Add meaningful names for the output
    )) %>% 
    arrange(desc(mean))
  
  # Tukey's test
  tukey.result <- TukeyHSD(anova.result)
  # print(tukey.result)
  
  # creating the compact letter display
  cld.result <- multcompLetters4(anova.result, tukey.result)
  # print(cld.result)
  
  # adding the compact letter display to the table with means and sd
  cld <- as.data.frame.list(cld.result$`Slaughter_Condition:Feeding_Regime`)
  Data_summary$Tukey <- cld$Letters
  # print(Data_summary)
  
  # Saving the excel file (Not able to save files)
  # write.csv(Data_summary, "data_ANOVA and Tukey_summary.csv", row.names = FALSE)
  
  ### Bar plot
  # x label and y label
  xlabel = "Feeding Regime"
  ylabel = ""
  
  # colored barplot according Tukey's test results
  plt <- ggplot(Data_summary, aes(x = Feeding_Regime, y = mean, group = Slaughter_Condition, col = Slaughter_Condition)) + 
    geom_point(aes(shape = Slaughter_Condition))+
    geom_line() +
    geom_errorbar(aes(ymin =  mean-sd, ymax = mean+sd), show.legend = FALSE) +
    # geom_bar(stat = "identity", width=0.8, alpha=0.8, show.legend = FALSE) +
    # scale_fill_brewer(palette = "BrBG") +
    # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.2) + 
    labs(x= xlabel, y= ylabel) +
    ggtitle(glue("{response_colname[i]}")) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_text(aes(label=Tukey), nudge_y = 0.1,size = 3) +
    scale_color_manual(values = c("NoStress" = "blue", "Stress" = "red"))
  
  # plt
  # saving the final figure
  ggsave(glue("plots/checking_interac/{i}_{response_colname[i]}_lineplot.png"), height = 8, width = 12, dpi = 1000)
  
  Data_summary$new <- paste(Data_summary$Slaughter_Condition, Data_summary$Feeding_Regime, sep = ".")
  plt2 <- ggplot(Data_summary, aes(new, mean, fill = Tukey)) + 
    geom_bar(stat = "identity", width=0.8, alpha=0.8, show.legend = FALSE) +
    scale_fill_brewer(palette = "BrBG") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.2) + 
    labs(x= xlabel, y= ylabel) +
    ggtitle(glue("{response_colname[i]}")) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_text(aes(label=Tukey), nudge_y = 0.1,size = 3)
  
  ggsave(glue("plots/checking_interac/{i}_{response_colname[i]}_barplot.png"), height = 8, width = 12, dpi = 1000)
  # dev.off()
}

# # ?ggsave
# 
# plt <- ggplot(Data_summary, aes(factor(), mean, fill = Tukey)) + 
#   geom_bar(stat = "identity", width=0.8, alpha=0.8, show.legend = FALSE) +
#   scale_fill_brewer(palette = "BrBG") +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.2) + 
#   labs(x= xlabel, y= ylabel) +
#   ggtitle(glue("{response_colname[i]}")) +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_text(aes(label=Tukey), nudge_y = 0.1,size = 3)



# ####
# library(readxl)
# upregulated_genes <- read_excel("data/07_01_01_up_and_down_regulated_genes.xlsx", sheet = "upregulated_genes")
# downregulated_gene <- read_excel("data/07_01_01_up_and_down_regulated_genes.xlsx", sheet = "downregulated_gene")
# 
# upregulated_genes$Accession[!upregulated_genes$Accession %in% response_colname]
