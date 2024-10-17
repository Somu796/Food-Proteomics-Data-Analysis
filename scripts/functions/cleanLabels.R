cleanLabels <- function(ColNames, Class, label_mappings) {
  NColNames <- c()
  
  for (i in 1:length(ColNames)) {
    new_label <- NA  # Default to NA if no match is found
    
    for (mapping in label_mappings[[Class]]) {
      if (grepl(mapping$pattern, ColNames[i])) {
        new_label <- mapping$new_label
        break
      }
    }
    
    NColNames <- append(NColNames, new_label)
  }
  
  return(NColNames)
}

# view(data.frame(new_lab = cleanLabels(ColNames, "Stress", label_mappings), old_lab = ColNames))



# CleanLabels(Normalized_Data$Label, "Stress")
# CleanLabels <- function(ColNames, Class) {
#   NColNames = c()
#   for (i in 1: length(ColNames)){
#     
#     if (grepl("Horn", Class)){
#       NColNames = append(NColNames, if (grepl("NC", ColNames[i])) "NoHorn" else "Horn") # For Horn
#       
#     }else if(grepl("Stress", Class)){
#       NColNames = append(NColNames, if (grepl("NoStress", ColNames[i])) "NoStress" else "Stress") # For Stress
#       
#     }else if(grepl("Rearing", Class)){
#       NColNames = append(NColNames, if (grepl("NM", ColNames[i])) "NoMixed" else "Mixed") # For Rearing
#       
#     }else{
#       NColNames = append(NColNames, NA)
#       
#     }
#   }
#   
#   return (NColNames)
# }





