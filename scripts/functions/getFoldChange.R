getFoldChange = function(AOVResultData, columns){
  
  FOldChangeData = AOVResultData %>%
    dplyr::select(all_of(columns)) %>% #Select only works on character variable
    filter(., !is.na(.[[Factor]]))
  
  FOldChangeData$MeanRatio = FOldChangeData[[meanFactorLevel.1]]/FOldChangeData[[meanFactorLevel.2]] #Calculating Mean ratio
  FOldChangeData$log2FC = log2(FOldChangeData$MeanRatio)
  
  return(FOldChangeData)
}

