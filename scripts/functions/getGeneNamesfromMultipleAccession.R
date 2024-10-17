# information = "accession,gene_primary,gene_synonym,organism_name,protein_name,sequence" #Information to be fetched from database 
# #Possible options: "accession,id,gene_names,gene_primary,gene_synonym,sequence,gene_oln,gene_orf,organism_name,organism_id,protein_name,xref_proteomes,lineage,virus_hosts"


getGeneNamesfromMultipleAccession <- function(information, gene_data) {
  
  UniprotNames = data.frame(matrix(nrow = 0, ncol = length(unlist(strsplit(information, ","))))) # Dummy data.frame to collect information from uniprot
  colnames(UniprotNames) = unlist(strsplit(information, ",")) # renaming the dummy data.frame colnames with the information header
  
  NA_GeneData = gene_data
  
  # Handling multiple accession number:
  no_accession = sum(grepl("Accession", colnames(gene_data)))
  
  for (n in 1:(no_accession-1)){
    print(n)
    UniprotNames_dummy = GetGeneUniProt(as.vector(NA_GeneData[[paste0("Accession.", n)]]), information) #Enriching Gene Information from uniprot using Accession Number
    colnames(UniprotNames_dummy) = str_to_title(unlist(strsplit(information, ","))) #Changing column to more readable names and capitalize first letter
    rownames(UniprotNames_dummy) = NULL #Rownames converted to NULL
    
    UniprotNames = rbind(UniprotNames, UniprotNames_dummy)
    
    ## If Accession.1 UniProt is NA (UniProt --> if NA)
    NA_geneprimary = c()
    
    NA_geneprimary = which(is.na(UniprotNames$Gene_primary)) # This will return entries where Accession.1 had NA, Accession.2, Accession.n
    
    ## Filtering gene_data according to Accession.1 UniProt information is NA (if NA --> Accession.n+1)
    NA_GeneData = c()
    NA_GeneData = gene_data[gene_data[[paste0("Accession.", n)]] %in% UniprotNames$Accession[NA_geneprimary],] 
    
    ## Filtering Accession numbers which are not NA in Acc.2 (Accession.n+1 --> if NOT NA)
    NA_GeneData  = NA_GeneData[!is.na(NA_GeneData[[paste0("Accession.", n+1)]]), ]
    
    # %>%
    #   filter(!is.na(Accession.2)) 
    
    if(nrow(NA_GeneData)==0){
      print("I entered")
      break
    }
    
    ## Dropping UniProt Entries (Accession.n+1 --> else drop entries in UniProt) 
    UniprotNames = UniprotNames %>%
      filter(Accession %in% setdiff(UniprotNames$Accession, NA_GeneData[[paste0("Accession.", n)]]))
    print(n)
  }
  return(UniprotNames)
}