app_match_mean_long <- function(data,clinic.dat) {
  
  prematched <- merge(data, clinic.dat, by = c("PID", "Date"), all.x = TRUE)
  
  unmatched <- prematched %>% filter(is.na(Diagnosis) == TRUE)
  
  list.id <- as.data.frame(unmatched$PID)
  
  list.date <- as.data.frame(unmatched$Date)
  
  mismatched <- bind_rows(lapply(1:nrow(list.id), function(x) match_clin(list.id[x,1], clinic.dat, data, list.date[x,1])))
  
  mismatched$date.diff <- as.numeric(mismatched$date.diff)
  
  mismatched <- as.data.frame(mismatched)
  
  prematched <- prematched %>% mutate(date.diff = 0) %>% filter(is.na(Diagnosis) == FALSE)
  
  prematched <- as.data.frame(prematched)
  
  print(nrow(mismatched))
  print(nrow(prematched))
  
  
  final <- bind_rows(prematched,mismatched) %>% dplyr::select(-date.diff)
  
  return(final)
  
}
