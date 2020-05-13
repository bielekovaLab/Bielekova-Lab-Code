find_vals <- function(data, chunk, cohort) {
  
  
  dat <- data %>% group_by(PID,Date) %>% filter(Groups == "90sec") %>% 
    mutate(tries_number = length(unique(Time))) %>% filter(tries_number == 2) %>%
    mutate(Trial_Num = ifelse(Time == min(Time), "Trial 1", "Trial 2")) %>%
    group_by(PID, Date, Trial_Num) %>% mutate(correct.count = ifelse(guessedid == correctid,
                                                                     1, 0)) %>% filter(Diagnosis == cohort)
  
  time.chunk <- dat %>% group_by(PID,Date,Trial_Num) %>% filter(Time.Vals/max(Time.Vals) <= chunk) %>%
    summarise(sum(correct.count ))
  
  casted.chunk <- cast(time.chunk, PID + Date ~ Trial_Num)
  
  return(casted.chunk)
  
}  
