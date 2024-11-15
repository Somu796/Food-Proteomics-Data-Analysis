# Description
# To Run: Run "01_cleaning_raw_data" scripts to run this.
# saving_data: yes
# data_name: "data/03_enriched_gene_name_information_from_uniprot_post_02_enriching_.RData", 
# 
# "data/05_merged_accession_enriched_gene_name_information.RData", 
# 
# "data/06_gene_data_post_02_01_with_multple_accession_entries.RData", 
# 


# 1. Libraries ------------------------------------------------------------
pacman::p_load(tidyverse, tidylog)
source("scripts/00_path_variables.R") #only run when you are not running, "01_cleaning_raw_data"

# library(UniprotR) # We have modified the source code because it has time limit on API call
source(countAccession)
source(getAccession)
source(getGeneNames)
source(getGeneNamesfromMultipleAccession)
source("scripts/functions/getMergedGeneNames.R")

# 2. Input for generating gene_data_01 from gene_data --------------------------------
# Input
# gene_data <- c("data/01_gene_data_post_01_cleaning.RData")
# Output
# save_gene_data_01_path <- c("data/06_gene_data_post_02_01_with_multple_accession_entries.RData")
# save_UniprotNames_path <- c("data/03_enriched_gene_name_information_from_uniprot_post_02_enriching_.RData")
# save_merged_Accession_UniprotNames_information_path <- c("data/05_merged_accession_enriched_gene_name_information.RData")

# 3. loop for generating gene_data_01 from gene_data ----------------------

load(gene_data) 

# 4. Preparing gene_data for UniProt Enrichment ---------------------------

## 4.1. There are multiple Accession number are there in few entries, so that has to be divided separated in multiple col names (modify) --------------


colName = "Accession"
delimiter = ";"

gene_data = gene_data %>%
  separate_wider_delim(Accession, ";", names = paste0("Accession.", 1:countAccession(gene_data, colName, delimiter, endsWith = FALSE)), too_few = "align_start")

# 4.2. Fetching the colnames have Accession in it for further cleaning of accession numbers ----------------------------
accession_colnames = colnames(gene_data)[grepl("Accession", colnames(gene_data))] 

# 4.3. Cleaning all the 9 accession number column with getAccession Function, eg: A5D7E1|A5D7E1_BOVIN -> A5D7E1; --------------------------------  
# source(getAccession)
for (i in accession_colnames){
  gene_data[[i]] = getAccession(gene_data, i)
}

# 5. Getting Gene Name both from Description and UniProt ---------------------------------
# source(getGeneNames)

## 5.1. Cleaning Gene Name from "Description"/"GeneName.Desc" ----------------------------------

colnames(gene_data)[which(colnames(gene_data) == "Description")] = "GeneName.Desc" # renaming column "Description" -> "GeneName.Desc"
gene_data$GeneName.Desc = GeteGeneNameDesc(gene_data, "GeneName.Desc")

## 5.2. Enrichment of Gene Name with UniProt Database ------------------------

# source(getGeneNamesfromMultipleAccession)

information = "accession,gene_primary,gene_synonym,organism_name,protein_name,sequence" #Information to be fetched from database
#Possible options: "accession,id,gene_names,gene_primary,gene_synonym,sequence,gene_oln,gene_orf,organism_name,organism_id,protein_name,xref_proteomes,lineage,virus_hosts"

UniprotNames <- getGeneNamesfromMultipleAccession(information, gene_data)

## 5.3. Enrichment of Gene Name with UniProt Database of merged_Accession  ------------------------------------------

### 5.3.1. enriching "merged" accessions with UniProt -------------------------
# if Protein_name is merged in uniprot data that means 

# source("scripts/functions/getMergedGeneNames.R")

merged_Accession_UniprotNames <-  UniprotNames %>% 
  filter(Protein_name == "merged") # filtering out only merged accession numbers

merged_Accession_UniprotNames_information = getMergedGeneNames(merged_Accession_UniprotNames$Accession) #collecting merged_accession -> new_accession

merged_Accession_UniprotNames_information_dummmy = GetGeneUniProt(as.vector(merged_Accession_UniprotNames_information[["new_Accession"]]), information)  #collecting details corresponding the merged accession number
colnames(merged_Accession_UniprotNames_information_dummmy) = str_to_title(unlist(strsplit(information, ","))) #Changing column to more readable names and capitalize first letter
rownames(merged_Accession_UniprotNames_information_dummmy) = NULL #Rownames converted to NULL

merged_Accession_UniprotNames_information <- merged_Accession_UniprotNames_information %>% 
  full_join(merged_Accession_UniprotNames_information_dummmy, by = join_by(new_Accession == Accession)) # making a join of [Accession, new_Accession] with UniProt information


### 5.3.2. Replacing merged merged_accession_uniprot with UniProtNam --------

UniprotNames <- UniprotNames %>% 
  filter(Protein_name != "merged") # removing the (Protein_names == merged) entries from UniprotNames

UniprotNames <- rbind(UniprotNames, merged_Accession_UniprotNames_information %>% 
                        select(-new_Accession)) #adding up the Uniprot information for merged entries to UniprotNames

# 6. Saving the data in .RData format for future use ----------------------------------------------------------
save(UniprotNames, file = save_UniprotNames_path)
save(merged_Accession_UniprotNames_information, file = save_merged_Accession_UniprotNames_information_path)

gene_data_01 = gene_data # this will make easy for 02_02_merging_gene data, I don't have to rerun the same analysis again
save(gene_data_01, file = save_gene_data_01_path)






