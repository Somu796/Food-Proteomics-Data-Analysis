getAccession <- function(Data, ColName){
  accessionNo = c()
  for (i in 1:nrow(Data)){
    accessionNo = append(accessionNo, sub("\\|.*", "", Data[i,ColName]))
  }
  if (sum(is.na(accessionNo)) == 0) print("No Na present in accession Numbers") else print("NA present in Accession Numbers")
  return(accessionNo)
}

# checking if any access number is na