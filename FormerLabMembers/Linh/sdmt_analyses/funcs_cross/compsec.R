compsec <- function(data) {
  
  time_score <- data %>% filter(Groups == "90sec") %>% group_by(PID, Date, Time) %>%
        mutate(count = ifelse(guessedid == correctid, 1, 0),
               psec_90 = Correct/(max(Time.Vals)/1000)) %>%
        filter(Time.Vals <= (75/90)*max(Time.Vals)) %>%
        group_by(PID, Diagnosis, Date, Time, psec_90) %>% 
        summarise(psec_75 = sum(count)/(max(Time.Vals)/1000)) %>%
        group_by(PID, Date) %>% mutate(`75 sec` = mean(psec_75),
                                             `90 sec` = mean(psec_90)) %>%
        distinct(`75 sec`, .keep_all = TRUE)
  
  ms <- time_score %>% filter(Diagnosis == "MS") %>% distinct(PID)
  hv <- time_score %>% filter(Diagnosis == "HV") %>% distinct(PID)
  
  p <- ggpaired(time_score, cond1 = "75 sec", cond2 = "90 sec",
           line.color = "gray", line.size = 0.4,
           title = "Points/Sec at 75 Seconds v. 90 Seconds") +
  border("black") + theme(axis.ticks.length=unit(-0.25, "cm"),
  axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  axis.title.x = element_text(size = 14, face = "bold"),
  axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  plot.title = element_text(hjust = 0.5)) + labs(y = "Points/Sec", x = "") +
  ylim(0.15,1.15)
  
  print(p)
  
  speed_range <- range(time_score$`90 sec`, time_score$`75 sec`)
  
  print(wilcox.test(time_score$`75 sec`, time_score$`90 sec`, paired = TRUE))
  print(paste("mean speed at 90s - at 75s is",c(mean(time_score$`90 sec`) - mean(time_score$`75 sec`))))
  print(paste("range of all speed at 75 and 90s is", speed_range[1], "to", speed_range[2]))
  print(paste("HV =", nrow(hv),"MS =", nrow(ms)))
  
}
