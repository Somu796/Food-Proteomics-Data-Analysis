# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: data/04_unused_accessions_for_enriching_gene_data.RData",
# "data/LT_data/04_unused_accessions_for_enriching_gene_data.RData"

# 1. Libraries ------------------------------------------------------------

pacman::p_load(tidyverse, glue, tidylog)

source("scripts/00_path_variables.R")


# 2. Merging gene_data with UniProt_data ----------------------------------

# Description
# To RUN: This is an independent script.
# saving_data: yes
# data_name: "data/04_unused_accessions_for_enriching_gene_data.RData",
#
# "data/07_enriched_gene_data_01_after_02_02_merging_gene_data_with_uniprot.RData",
# 

# 1. Libraries ------------------------------------------------------------

library(tidyverse)
library(glue)
library(tidylog)

source("scripts/00_path_variables.R")

#input
# gene_data_01 <- c("data/06_gene_data_post_02_01_with_multple_accession_entries.RData")
# UniprotNames <- c("data/03_enriched_gene_name_information_from_uniprot_post_02_enriching_.RData")

#output
# save_unused_accession_path <- c("data/04_unused_accessions_for_enriching_gene_data.RData")
# save_enriched_gene_data_01_path <- c("data/07_enriched_gene_data_01_after_02_02_merging_gene_data_with_uniprot.RData")


# 2. Merging gene_data with UniProt_data ----------------------------------

## 2.1. Loading Data -------------------------------------------------------
load(gene_data_01)
load(UniprotNames) ## Data is saved in Uniprot

## 2.4. wide_to_longer: Converting Accession.1,..., Accession.9 [863*89]-> Accession of longer data_inputs [1121*81]-----------------

gene_data_01 = gene_data_01  %>%
  # converting all the accession as individual entries
  pivot_longer(
    cols = starts_with("Accession."), 
    names_to = "Accession_no", #creating extra column "Accession_no" to saving number Accession."1, 2, 3, 4"
    names_prefix = "Accession.",#to remove "Accession." from Accession.1, Accession.2,... 
    values_to = "Accession",
    values_drop_na = TRUE
  ) %>%
  relocate(c("Accession", "Accession_no"), .before = 1)

## 2.5. table for multiple_accession: how many entries has how many accession ----------------------------

accession = c()
entries_count = c()
multiple_entries = c()
j = 0 # to ignore sum(Accession.9) - sum(accession.10)

for (i in max(gene_data_01$Accession_no):min(gene_data_01$Accession_no)){
  accession = append(accession, i) 
  entries_count = ifelse(j>0, 
                         glue("Only {sum(gene_data_01$Accession_no == i) - sum(gene_data_01$Accession_no == (i+1))} entries has {i} accessions."),
                         glue("Only {sum(gene_data_01$Accession_no == i)} entries has {i} accessions."))
  
  multiple_entries = append(multiple_entries, entries_count)
  j = j+1
}


## 2.6. Removing extra column and saving unused Accessions in a list  --------------

### 2.6.1. Removing extra column "Accession_no"----------
gene_data_01 = gene_data_01 %>%
  select(-"Accession_no")

### 2.6.2. saving unused accession numbers -------------
unused_accession = data.frame(unused_accession = setdiff(gene_data_01$Accession, UniprotNames$Accession)) #saving unused accessions for enriching gene data
save(unused_accession, file = save_unused_accession_path)

## 2.7. Merging gene_data_01 with UniprotNames --------------------------------
enriched_gene_data =  UniprotNames %>%
  left_join(gene_data_01, by= join_by(Accession)) %>% 
  relocate(c("GeneName.Desc"), .before = "Organism_name") # left join of uniprot to gene_data_01

## 2.8. saving the enriched_data (There are multiple possibilities) -------

### 2.8.1. Final Decision: Focus only on Gene.primary and Protein_name != "deleted, otherwise others might have been modified, with acccession number, I can support better
enriched_gene_data_01 = enriched_gene_data[!(is.na(enriched_gene_data$Gene_primary) & enriched_gene_data$Protein_name == "deleted"), ] # If both description and Uniprot has NA 
save(enriched_gene_data_01, file = save_enriched_gene_data_01_path)


# 3. output to show ---------------------
multiple_entries
