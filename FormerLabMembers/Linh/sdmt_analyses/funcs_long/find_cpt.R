find_cpt <- function(data, patient, c_input, b0_input, b1_input) {
  
  pat <- patient
  y <- data %>% filter(PID==pat) %>% .$y
  x <- data %>% filter(PID==pat) %>% .$x
  
  mod <- nls(y~fx(x,c,b0,b1), start=c("c"= c_input,"b0"= b0_input,"b1"= b1_input),
             control=nls.control(warnOnly = TRUE,minFactor = 1e-20,maxiter=1000))
  
  out <- summary(mod)
  
  
  fin <- data.frame(rep(pat, length(y)), x, y, predict(mod),
             rep(round(out$coefficients[1]), length(y)))
  
  colnames(fin) <- c("PID", "x", "y", "pred_y", "cpt")
  
  return(fin)
  
}




