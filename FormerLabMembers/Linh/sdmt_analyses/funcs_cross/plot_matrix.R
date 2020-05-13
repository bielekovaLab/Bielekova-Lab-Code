plot_matrix <- function(input, test) {
  
  dat <- fig5
  
  dat.elim <- dat %>% dplyr::select(-PID,-Date,-Diagnosis,-Gender,-PASAT,-App,-Written,
                                    -`Upper Strength`)
  
  dat.elim <- as.data.frame(dat.elim)
  
  if(test == "PASAT") {
    y <- dat$PASAT
  }
  
  else if(test == "App") {
    y <- dat$App
  }
  else{
    y <- dat$Written
  }
  
  par(mfrow = c(4,3), mar = c(5,5,1,2))
  
  plot_all <- function(x) { 
    
  plot(y~dat.elim[,x], xlab = "",
       ylab = "")
    
  mtext(colnames(dat.elim)[x], side = 1, line = 3)  
  mtext(test, side = 2, line = 3)
  
  abline(lm(y~dat.elim[,x]), col = "red")
  
  }
  
  sapply(1:ncol(dat.elim), function(x) plot_all(x))
  
  
}
  