outlier.adj <- function(input,option) {
  
  thresh.list <- mix_ba1v2(input, plot_option = FALSE)
  
  thresh.lo <- thresh.list[1]
  thresh.mean <- thresh.list[2]
  thresh.hi <- thresh.list[3]
  
  
    data <- input %>% group_by(PID,Date,Time) %>% 
      filter(Time.Vals == min(Time.Vals)) %>%
      ungroup()
    
    outlier.adj <- data %>% group_by(PID) %>% group_by(PID,Date) %>% mutate(Correct = score_p_sec*90) %>%
      mutate(diff = Correct[which(Time == max(Time))] - Correct[which(Time == min(Time))]) %>%
      mutate(Correct.adj = ifelse(diff > thresh.hi & Correct == min(Correct), 
                                     Correct + (diff - thresh.hi),
                           ifelse(diff < thresh.lo & Correct == min(Correct), 
                                         Correct + abs(diff - thresh.lo),
                                         Correct))) 
    
    test <- outlier.adj %>% group_by(PID, Date) %>% filter(n() == 2) %>% 
            summarise(badiff = Correct.adj[which(Time == max(Time))] - Correct.adj[which(Time == min(Time))],
                      mean = (Correct.adj[which(Time == max(Time))]+Correct.adj[which(Time == min(Time))])/2) 
    
   
    outliers.up <- outlier.adj %>% group_by(PID, Date) %>% 
                 filter(diff > thresh.hi) %>% group_by(PID, Date) %>% filter(n() == 2) %>% 
                summarise(badiff = Correct.adj[which(Time == max(Time))] - Correct.adj[which(Time == min(Time))],
                mean = (Correct.adj[which(Time == max(Time))]+Correct.adj[which(Time == min(Time))])/2) 
    
    outliers.lo <- outlier.adj %>% group_by(PID, Date) %>% 
      filter(diff < thresh.lo) %>% group_by(PID, Date) %>% filter(n() == 2) %>% 
      summarise(badiff = Correct.adj[which(Time == max(Time))] - Correct.adj[which(Time == min(Time))],
                mean = (Correct.adj[which(Time == max(Time))]+Correct.adj[which(Time == min(Time))])/2) 
    
    if(option == TRUE) { 
    plot(badiff~mean, data = test, ylim = c(-25,25), xlim = c(10,110), tck = 0.02)
    abline(h = thresh.lo, lty = 2, lwd = 2)
    abline(h = thresh.hi, lty = 2, lwd = 2)
    abline(h = 0, lwd = 2)
    abline(h = thresh.mean, lty = 2, lwd = 2, col = "red")
    points(outliers.up$mean, outliers.up$badiff, col = "blue")
    points(outliers.lo$mean, outliers.lo$badiff, col = "chocolate3")
    
    }
    
    if(option == FALSE){
    
    fin <- outlier.adj %>% mutate(score_p_sec = Correct.adj/90) %>% 
           dplyr::select(-Correct.adj, -diff)
    
    return(fin)
      
    }
  
  
  
    
}
