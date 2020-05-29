# these are functions to compare immune deconvolution resuilts
# the immune deconv input data frame must have cell types as row names and patient IDs as column names
# this MUST be looped into a a tumor_types variable
# tumor_types <- c('BRCA','UCEC','OV','THCA','LGG','SARC','COAD')


deconvCompareStat <- function(deconv_matrix, surv_data = survival_data, statistic){
  
  # compare extreme quartiles with either wilcox or t test
  
  age_subset <- survival_data %>%
    dplyr::filter(type %in% tumor_types[i],
                  Quartiles %in% c(1, 4)) %>%
    dplyr::select(bcr_patient_barcode, Quartiles) %>%
    dplyr::mutate(Quartiles = ifelse(Quartiles == 1, "Young", "Old"))
  age_subset <- age_subset[!is.na(age_subset$Quartiles), ]
  
  cases <- intersect(colnames(deconv_matrix), age_subset$bcr_patient_barcode)
  
  temp_xcell <- deconv_matrix %>% 
    dplyr::select(cases)
  temp_xcell <- temp_xcell[, sort(colnames(temp_xcell))]
  
  age_subset <- age_subset %>%
    dplyr::filter(bcr_patient_barcode %in% cases) %>%
    dplyr::arrange(bcr_patient_barcode)
  
  if(!all(colnames(temp_xcell) == age_subset$bcr_patient_barcode)){print("ERROR: ", tumor_types[i])}
  
  tmp <- merge(t(temp_xcell), age_subset, by.x = 0, by.y = "bcr_patient_barcode")
  tmp <- data.table::melt(tmp)
  
  if (statistic == "wilcox"){
    wilcox_result <- tmp %>%
      dplyr::group_by(variable) %>%
      rstatix::wilcox_test(value ~ Quartiles, exact = F, detailed = T)
    
    young_med <- tmp %>%
      dplyr::group_by(variable, Quartiles) %>%
      dplyr::filter(Quartiles == "Young") %>%
      dplyr::summarise(Young_Median = median(value))
    old_med <- tmp %>%
      dplyr::group_by(variable, Quartiles) %>%
      dplyr::filter(Quartiles == "Old") %>%
      dplyr::summarise(Old_Median = median(value))
    med_scores <- data.frame(Cell = young_med$variable,
                             Young_Median = young_med$Young_Median,
                             Old_Median = old_med$Old_Median) %>%
      dplyr::mutate(Higher_Med = ifelse(Young_Median > Old_Median, "Young",
                                        ifelse(Young_Median == Old_Median, "No Difference", "Old")))
    results <- merge(wilcox_result, med_scores, by.x = "variable", by.y = "Cell")
  }
  
  if (statistic == "t.test"){
    results<- tmp %>%
      dplyr::group_by(variable) %>%
      rstatix::t_test(value ~ Quartiles, detailed = T)
  }
  
  
  return(results)
}


deconvCorrelateStat <- function(deconv_matrix,
                                surv_data = survival_data,
                                statistic = "pearson", comparision = "quartiles"){
  
  if(comparision == "quartiles"){
    age_subset <- survival_data %>%
      dplyr::filter(type %in% tumor_types[i],
                    Quartiles %in% c(1, 4)) %>%
      dplyr::select(bcr_patient_barcode, type, Quartiles) %>%
      dplyr::mutate(Quartiles = ifelse(Quartiles == 1, 0, 1))
    
    cases <- intersect(colnames(deconv_matrix), age_subset$bcr_patient_barcode)
    
    temp_xcell <- deconv_matrix %>% 
      dplyr::select(cases)
    temp_xcell <- temp_xcell[, sort(colnames(temp_xcell))]
    
    age_subset <- age_subset %>%
      dplyr::filter(bcr_patient_barcode %in% cases) %>%
      dplyr::arrange(bcr_patient_barcode)
    
    if(!all(colnames(temp_xcell) == age_subset$bcr_patient_barcode)){print("ERROR")}
    
    results <- apply(temp_xcell, 1,
                     function(x) cor.test(x, age_subset$Quartiles,
                                          method = statistic,
                                          exact = FALSE))
  }
  
  if(comparision == "continuous"){
    age_subset <- survival_data %>%
      dplyr::filter(type %in% tumor_types[i]) %>%
      dplyr::select(bcr_patient_barcode, age_at_initial_pathologic_diagnosis)
    
    cases <- intersect(colnames(deconv_matrix), age_subset$bcr_patient_barcode)
    
    temp_xcell <- deconv_matrix %>% 
      dplyr::select(cases)
    temp_xcell <- temp_xcell[, sort(colnames(temp_xcell))]
    
    age_subset <- age_subset %>%
      dplyr::filter(bcr_patient_barcode %in% cases) %>%
      dplyr::arrange(bcr_patient_barcode)
    
    if(!all(colnames(temp_xcell) == age_subset$bcr_patient_barcode)){print("ERROR")}
    
    results <- apply(temp_xcell, 1,
                     function(x) cor.test(x, age_subset$age_at_initial_pathologic_diagnosis,
                                          method = statistic,
                                          exact = FALSE))
  }
  
  return(results)
}
