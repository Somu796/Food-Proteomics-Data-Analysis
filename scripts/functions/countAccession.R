# objective: counting number of accession with help of the delimiter
## gene_data : data frame having the Accession column, e.g., "gene_data"
## colName : accession column name, e.g., "Accession"
## delimiter : what is the delimiter separating the accession numbers, e.g., ";" with special character add // at the beginning
## endsWith : does the accessions end with delimiter, default is FALSE, in case of TRUE, e.g., "A0A3Q1M2L8|A0A3Q1M2L8_BOVIN;"

pacman::p_load(glue)

countAccession <-  function(gene_data, colName, delimiter, endsWith = FALSE){
  max_count <- 0
  for (i in 1:length(gene_data[[colName]])){
    # delimiter <- glue("\\{delimiter}")
    matches <- gregexpr(delimiter, gene_data[[colName]][i])[[1]]
    count <- ifelse(matches[1] == -1, 0, length(matches))
    
    if (endsWith){
      # print(count)
      if (max_count < count){
        max_count <- count
        index <- i
      }
    }else{
      count <- count + 1
      # print(count)
      if (max_count < count){
        max_count <- count
        index <- i
      }
    } 
    
  }
  return(max_count)
}


# # load(gene_data_LT)
# 
# colName <- "Accession" 
# delimiter <- ";"
# a <- countAccession(gene_data, colName, delimiter, endsWith = FALSE)

