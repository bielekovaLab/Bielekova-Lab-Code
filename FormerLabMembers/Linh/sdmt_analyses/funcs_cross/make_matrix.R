make_matrix <- function(input, plot_option) {
  
  dat <- fig5
  
  hv <- dat %>% filter(Diagnosis == "HV") %>% distinct()
  ms <- dat %>% filter(Diagnosis == "MS") %>% distinct()
  print(paste("HV =", nrow(hv),"MS =", nrow(ms)))
  
  if(plot_option == TRUE) { 
  
  
  dat.elim <- dat %>% dplyr::select(-PID,-Date,-Diagnosis,-Gender,-`Upper Strength`)
  
  corr <- round(cor(dat.elim, use = "complete.obs", method = "spearman"),2)
  corr.adj <- reshape2::melt(corr)
  
  p.mat <- cor_pmat(dat.elim, method = "spearman")
  p.melt <- reshape2::melt(p.mat)

  fin <- merge(corr.adj, p.melt, by = c("Var1", "Var2")) %>%
    filter(Var1 != "App" &
             Var1 != "Written" &
             Var1 != "PASAT") %>%
    filter(Var2 == "PASAT" |
             Var2 == "Written" |
             Var2 == "App")

  label.assign <- fin %>% mutate(stars = ifelse(value.y < 0.01 & value.y >= 0.001, "*",
                                                ifelse(value.y < 0.001 & value.y >= 0.0001, "**",
                                                       ifelse(value.y < 0.0001, "***", "")))) %>%
    mutate(val = ifelse(value.x == 1, "",
                        paste(stars, formatC(round(value.x,2),
                                             format = 'f', digits = 2), sep="\n"))) %>% 
            dplyr::rename(`Correlation` = value.x) 
  
  label.assign$Var1 <- factor(label.assign$Var1,
                                     levels = c("Brain Volume",
                                                "Lesions Volume",
                                                "Cognitive Functions", 
                                                "DH Cerebellar",
                                                "DH Sensory", "Brainstem/LC Nerves",
                                                "Eyes", "Eye Movement","EDSS","tap", "Age"))
  
  q <- ggplot(data = label.assign, aes(x=Var1, y=Var2, fill = Correlation, label= Correlation)) 
  
  q + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank()) + geom_tile(colour = "gray",
                                                                          size = 0.25) +
    scale_fill_gradient2(low = "#0072B2",high = "#D55E00"
                         ,mid = "white") + geom_text(aes(label = label.assign$val),
                                                     hjust= 0.5, vjust = 0.5, size = 5)
  }
  
  else if(plot_option == FALSE) { 
  
  output <- paste("There are ", length(unique(subset(dat, Diagnosis == "HV")$PID)), "HV and",
                  length(unique(subset(dat,Diagnosis == "MS"))$PID), "MS")
  
  print(output)
  
  return(dat)
  }
  
  else (NULL)
}
