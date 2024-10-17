# It has two functions in it, 1. Get Gene name from description and Get Gene Name from UniProt

# 1. Get Gene Names from Description --------------------------------------
library(stringr)

GeteGeneNameDesc = function(Data, ColName){
  pattern <- "GN=([^[:space:]]+)" # Define a regular expression pattern to match the gene name
  gene_list = c() # Extract the gene name using str_extract
  for (i in 1:nrow(Data)){
    gene_name <- str_extract(Data[i,ColName], pattern) %>%
      {if (!is.na(.)) str_replace(., "GN=", "") else NA}
    gene_list = append(gene_list, gene_name)
  }
  if (sum(is.na(gene_list)) == 0) print("No Na present in gene list") else print("NA present in gene list")
  return(gene_list)
}


# 2. Get Gene Names from UniProt ------------------------------------------
library(httr)
GetGeneUniProt <- function(ProteinAccList,  ColName, directorypath = NULL)
{
  
  message("Please wait we are processing your accessions ...")
  pb <- progress::progress_bar$new(total = length(ProteinAccList))
  
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  #ColName = "accession,id,gene_names,gene_primary,gene_synonym,gene_oln,gene_orf,organism_name,organism_id,protein_name,xref_proteomes,lineage,virus_hosts"
  
  ProteinInfoParsed_total = data.frame()
  for (ProteinAcc in ProteinAccList)
  {
    #to see if Request == 200 or not
    Request <- tryCatch(
      {
        GET(paste0(baseUrl , ProteinAcc,"&format=tsv"))
      },error = function(cond)
      {
        message("Internet connection problem occurs and the function will return the original error")
        message(cond)
      }
    ) 
    #this link return information in tab formate (format = tab)
    ProteinName_url <- paste0(ProteinAcc,"&format=tsv&fields=",ColName)
    RequestUrl <- paste0(baseUrl , ProteinName_url)
    RequestUrl <- URLencode(RequestUrl)
    if (length(Request) == 0)
    {
      message("Internet connection problem occurs")
      return()
    }
    if (Request$status_code == 200){
      # parse the information in DataFrame
      ProteinDataTable <- tryCatch(read.csv(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
      if (!is.null(ProteinDataTable))
      {
        ProteinDataTable <- ProteinDataTable[1,]
        ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
        # add Dataframes together if more than one accession
        ProteinInfoParsed_total <- rbind(ProteinInfoParsed_total, ProteinInfoParsed)
      }
      
    }else {
      HandleBadRequests(Request$status_code)
    }
    pb$tick()
    
  }
  if(!is.null(directorypath))
  {
    write.csv(ProteinInfoParsed_total , paste0(directorypath ,"/" , "Names & Taxa Information.csv"))
  }
  return(ProteinInfoParsed_total)
}