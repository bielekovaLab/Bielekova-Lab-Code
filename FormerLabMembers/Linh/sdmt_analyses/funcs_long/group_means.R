group_means <- function(data, pt) {
  
  # filter for the patient's data from all longitudinal data
  patient <- data %>% filter(PID == pt)
  
  # calculate the mean scores for grouping of 1 sitting
  y <- rollapply(patient$score, 1, mean, by = 1, align = "left")
  x <- seq(1,length(y), 1)
  pid <- rep(pt, length(y))
  grouping <- rep(1, length(y))
  mean_1 <- data.frame(grouping,pid,x,y)
  
  # calculate the mean scores for grouping of 2 sittings
  y <- rollapply(patient$score, 2, mean, by = 2, align = "left")
  x <- seq(1,length(y), 1)
  pid <- rep(pt, length(y))
  grouping <- rep(2, length(y))
  mean_2 <- data.frame(grouping,pid,x,y)
  
  # calculate the mean scores for grouping of 3 sittings
  y <- rollapply(patient$score, 3, mean, by = 3, align = "left")
  x <- seq(1,length(y), 1)
  pid <- rep(pt, length(y))
  grouping <- rep(3, length(y))
  mean_3 <- data.frame(grouping,pid,x,y)
  
  # calculate the mean scores for grouping of 4 sittings
  y <- rollapply(patient$score, 4, mean, by = 4, align = "left")
  x <- seq(1,length(y), 1)
  pid <- rep(pt, length(y))
  grouping <- rep(4, length(y))
  mean_4 <- data.frame(grouping,pid,x,y)
  
  fin <- bind_rows(mean_1, mean_2,mean_3,mean_4)
  
  return(fin)
}

