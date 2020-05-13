app_match_mean <- function(data,clinic.dat) {
  
  
  prematched <- merge(data, clinic.dat, by = c("PID", "Date"))
  
  '%ni%' <- Negate('%in%')
  
  index <- which(data$PID %ni% prematched$PID)
  
  list.id <- as.data.frame(data$PID[index])
  
  list.date <- as.data.frame(data$Date[index])
  
  mismatched <- bind_rows(lapply(1:nrow(list.id), function(x) match_clin(list.id[x,1], clinic.dat, data, list.date[x,1])))

  mismatched$date.diff <- as.numeric(mismatched$date.diff)
  
  mismatched <- as.data.frame(mismatched)
  
  prematched <- prematched %>% mutate(date.diff = 0)
  
  prematched <- as.data.frame(prematched)
  
  final <- bind_rows(prematched,mismatched) %>% dplyr::select(-date.diff)
  
  return(final)
  
}
