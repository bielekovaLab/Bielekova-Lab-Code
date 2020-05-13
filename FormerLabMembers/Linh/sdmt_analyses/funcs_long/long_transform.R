long.transform <- function(input) {
  
  dat <- input
  
  dat <- dat %>% 
    dplyr::select(-Time) %>% 
    group_by(PID,Date,Diagnosis) %>% 
    summarise(y = mean(Correct)*90) %>% 
    ungroup() %>% 
    group_by(PID) %>% 
    mutate(x = rank(Date)) %>% 
    ungroup()
  
  return(dat)
}
