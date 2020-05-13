match_clin <- function(person,b_matched,app_dat, date) {
  
  patient <- app_dat %>% filter(PID == person & Date == date)
  
  clin <- b_matched %>% filter(PID == toString(person)) %>%
          mutate(date.diff = abs(Date - patient$Date[1])) %>%
          filter(date.diff == min(date.diff) & date.diff <= 3)
  
  if(nrow(clin) > 0) {
    
  clin$Date <- patient$Date[1]
  
  matched <- merge(patient, clin, by = c("PID", "Date"))
  
  return(matched)
  
  }
  else{
    NULL
  }
}
