multi_scale_match <- function(input,scale) { 
  
  cross.avg <- input %>% group_by(PID, Diagnosis,Date,Age,Gender) %>%
               summarise(score_p_sec = mean(score_p_sec))
  
  
  w.paper <- app_match_mean(cross.avg,paper) 
  w.pasat <- app_match_mean(cross.avg, pasat)
  w.tap <- app_match_mean(cross.avg, tap)
  
  if (scale == "mri") {
    df.list <- list(w.paper,w.pasat)
    merge.1 <- Reduce(function(x, y) merge(x, y, all = TRUE), df.list, accumulate=FALSE)
    w.scale <- app_match_mean(merge.1,mri)
  }
  
  else if (scale == "neurex") { 
    df.list <- list(w.paper,w.pasat)
    merge.1 <- Reduce(function(x, y) merge(x, y, all = TRUE), df.list, accumulate=FALSE)
    w.scale <- app_match_mean(merge.1,neurex) 
  } 
  
  else {
    w.scale <- app_match_mean(cross.avg, neurex_edss)
    
  }
  
  
  if(scale == "vis") {
    
    df.list <- list(w.paper,w.pasat)
    merge.1 <- Reduce(function(x, y) merge(x, y, all = TRUE), df.list, accumulate=FALSE)
    test <- app_match_vis(merge.1,vis) 
    
  }
  
  else { 
    df.list <- list(w.paper,w.pasat,w.scale,w.tap)
    
    test <- Reduce(function(x, y) merge(x, y, all = TRUE), df.list, accumulate=FALSE)
    
  }
  
  return(test)
  
}
