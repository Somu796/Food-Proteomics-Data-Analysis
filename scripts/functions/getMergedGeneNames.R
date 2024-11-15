# if (!require("httr")) install.packages("httr")
# if (!require("jsonlite")) install.packages("jsonlite")


library(httr)
library(jsonlite)

getMergedGeneNames <- function(Accession) {
  new_accession_list <- c()
  base_url <- "https://www.uniprot.org/uniprot/"
  query_url <- paste0(base_url, Accession, ".json")
  for (i in query_url) {
    response <- GET(i)
    if (status_code(response) != 200) {
      data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      if (!is.null(data$path)) {
        new_accession = sub(".*/", "", data$path)
        
      } else {
        new_accession = NA
      }
    } else {
      new_accession = NA
    }
    new_accession_list = append(new_accession_list, new_accession)
  }
  merged_accession = data.frame(Accession = Accession, new_Accession = new_accession_list)
  return(merged_accession)
}


# # Examples
# # Extract unique accession numbers
# Accession = c("F1MRX5", "F1N0W6", "F1MCF8")
# 
# # Get updated accession numbers and protein names
# 
# gerMergedGeneNames(Accession)
