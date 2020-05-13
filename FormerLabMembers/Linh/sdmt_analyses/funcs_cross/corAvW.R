corAvW <- function(data, type) {

matched <- fig4B

bamean <- baAvW(cross, option = FALSE)

adj <- matched %>% mutate(score_p_sec = score_p_sec + bamean)

hv <- adj %>% filter(Diagnosis == "HV") %>% distinct(PID)
ms <- adj %>% filter(Diagnosis == "MS") %>% distinct(PID)

print(paste("HV =", nrow(hv),"MS =", nrow(ms)))

lims <- range(adj$Paper_Correct, adj$score_p_sec, matched$Paper_Correct, matched$score_p_sec)

par(pty = 's')

if(type == "orig") { 
lims <- range(matched$Paper_Correct, matched$score_p_sec)

  plot(Paper_Correct~score_p_sec, data = matched,
       xlab = "App SDMT Score", 
       ylab = "Written SDMT Score",
       main = "App v. Written Score", col = "black", tck = 0.02,
       xlim = lims, ylim = lims)
  
  abline(0,1, lty = 1)
  abline(lm(Paper_Correct~score_p_sec, data = matched), col = "chocolate3", 
         lty = 2, lwd = 2)
  
  CCC <- CCC(matched$Paper_Correct, matched$score_p_sec)
  
  print(CCC$`rho.c`)
  print(summary(lm(Paper_Correct~score_p_sec, data = matched)))
}

if(type == "adj") {
  
  plot(Paper_Correct~score_p_sec, data = adj,
       xlab = "App SDMT Score", 
       ylab = "Written SDMT Score",
       main = "App v. Written Score", col = "black", tck = 0.02,
       xlim = lims, ylim = lims)
  
  abline(0,1, lty = 1)
  abline(lm(Paper_Correct~score_p_sec, data = adj), col = "chocolate3", 
         lty = 2, lwd = 2)
  
  CCC <- CCC(adj$Paper_Correct, adj$score_p_sec)
  
  print(CCC$`rho.c`)
  print(summary(lm(Paper_Correct~score_p_sec, data = adj)))
  
}
  print(nrow(matched))
  print(head(matched))
}
