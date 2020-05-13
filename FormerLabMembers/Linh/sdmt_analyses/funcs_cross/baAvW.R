baAvW <- function(data,option) { 

  matched <- fig4C
  
  orig.ba <- matched %>% mutate(badiff = Paper_Correct - score_p_sec*90,
                                bamean = (Paper_Correct + score_p_sec*90)/2)
  
  if(option == TRUE) { 
    
  par(pty = 's')
    
  plot(badiff~bamean, data = orig.ba, ylim = c(-15,30), tck = 0.02)
  
  abline(h = 0, lwd = 2)
  abline(h = mean(orig.ba$badiff), lty = 2,lwd=2, col = "red")
  abline(h = mean(orig.ba$badiff) + 1.96*sd(orig.ba$badiff), lty = 2,lwd=2)
  abline(h = mean(orig.ba$badiff) - 1.96*sd(orig.ba$badiff), lty = 2,lwd=2)
  
  print(mean(orig.ba$badiff) - 1.96*sd(orig.ba$badiff))
  print(mean(orig.ba$badiff))
  print(mean(orig.ba$badiff) + 1.96*sd(orig.ba$badiff))
  
  }
  
  if(option == FALSE) { 
 
  return(mean(orig.ba$badiff))
  }

}
