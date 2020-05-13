mix_ba_diag <- function(input) { 
  
  dat <- input %>% group_by(PID, Date, Time) %>% filter(Time.Vals == min(Time.Vals)) %>%
    ungroup() %>%
    group_by(PID, Date) %>% filter(n() == 2) %>%
    ungroup() %>% group_by(PID, Date, Diagnosis) %>% 
    summarise(mean_score = mean(score_p_sec)*90, 
              diff_score = (((score_p_sec[which(Time == max(Time))]))-
                              ((score_p_sec[which(Time == min(Time))])))*90) %>% 
    ungroup() %>% group_by(PID) %>% mutate(sitting = rank(Date) - 1) %>% ungroup()
  
  res <- lme(diff_score ~ sitting , random=~1|PID,
             correlation=corCompSymm(form=~1|PID), data = dat)
  
  totalsd<-sqrt(as.numeric(VarCorr(res)[1,1])+as.numeric(VarCorr(res)[2,1]))
  
  res2<-lme(diff_score~1,random=~1|PID,
            correlation=corCompSymm(form=~1|PID),data=dat)
  
  p.res <- plot(res)
  q.res <- qqnorm(res)
  
  p.res2 <- plot(res2)
  q.res2 <- qqnorm(res2)
  
  grid.arrange(p.res,q.res,p.res2,q.res2, ncol = 2, nrow = 2)
  
}
