all_long <- function(input,type) {
  
  if(type == "no age"){ 
  more20 <- input %>% group_by(PID, Date, Time) %>% filter(time.point == min(time.point)) %>% ungroup() %>% 
            group_by(PID) %>% mutate(number_o_tries = length(unique(Date))) %>% 
            filter(length(unique(Date)) >= 20) %>% ungroup() %>%
            dplyr::select(PID,Date,Time,score_p_sec,Diagnosis) %>% dplyr:: rename(Correct = score_p_sec)
  
  cat("There are", length(unique(subset(more20, Diagnosis == "HV")$PID)), "HV and", 
      length(unique(subset(more20, Diagnosis == "MS")$PID)), "MS")
  }
  
  else{
    more20 <- input %>% group_by(PID, Date, Time) %>% filter(time.point == min(time.point)) %>% ungroup() %>% 
      group_by(PID) %>% mutate(number_o_tries = length(unique(Date))) %>% 
      filter(length(unique(Date)) >= 20) %>% ungroup() %>%
      dplyr::select(PID,Date,Time,score_p_sec, Age,Diagnosis) %>% dplyr:: rename(Correct = score_p_sec)
  }
  
  return(more20)
}

