dat_part <- function(input, output_type, part_size) {

  #for analyzing with mri data (rm ACW083 because it's a technical issue)
  
  thing1 <- multi_scale_match(input,"mri") %>% dplyr::rename("Written" = Paper_Correct, "App" = score_p_sec) %>%
    filter(PID != "ACW083") %>% mutate(App = App*90)
  
  thing3 <- multi_scale_match(cross,"neurex") %>%
    dplyr::rename("Written" = Paper_Correct,
                  "App" = score_p_sec) %>% filter(PID != "ACW083") %>% mutate(App = App*90)
                  
  
  thing <- merge(thing1,thing3, by = c("PID","Date","Age","Diagnosis","PASAT","Written","App","Gender", "tap")) %>%
    na.omit() %>% group_by(Gender) %>% mutate(age_group = quant_bin(Age),
                                   app_group = quant_bin(App),
                                   writ_group = quant_bin(Written),
                                   pas_group = quant_bin(PASAT)) %>%
    ungroup() %>% mutate(sample_group = interaction(app_group,age_group, writ_group,pas_group)) %>%
    group_by(sample_group) %>%
    dplyr::select(-Volume_CerebellumGM, - Volume_CerebellumWM, -Volume_calculated_NAWMFr,-`DH Sensory`)
  
  set.seed(123)
  
  train <- sample_frac(thing, size = part_size) %>% ungroup %>% dplyr::select(-Date,-Gender,
                                                                  -age_group,-app_group,-sample_group,
                                                                  -writ_group,-pas_group)
  
  '%ni%' <- Negate('%in%')
  
  pid.set <- thing$PID[which(thing$PID %ni% train$PID)]
  
  test <- thing %>% filter(PID %in% pid.set) %>% ungroup %>% dplyr::select(-PID,-Date,-Gender,
                                                                    -age_group,-app_group,-sample_group,
                                                                    -writ_group,-pas_group)
  train <- train %>% dplyr::select(-PID)
  
  if (output_type == "train") {
    return(train)
  }
  
  else (return(test))
}
