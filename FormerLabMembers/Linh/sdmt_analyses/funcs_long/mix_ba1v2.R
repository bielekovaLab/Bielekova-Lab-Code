mix_ba1v2 <- function(input, plot_option) { 
  
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


mean<-summary(res2)$tTable[1,1]

# 95% Limits of agreement
low<-mean-1.96*totalsd
upper<-mean+1.96*totalsd

if(plot_option == TRUE) {

plot(diff_score~mean_score, data = dat, tck = 0.02, ylim = c(-25,25), xlim = c(10,110))

abline(h = 0, lwd = 2)
abline(h = mean, lty = 2, lwd = 2, col = "red")
abline(h = low, lty = 2, lwd = 2)
abline(h = upper, lty = 2, lwd = 2)

# get points that fall outside of the 95% limits of agreement

outliers.up <- dat %>% filter(diff_score > upper)
outliers.lo <- dat %>% filter(diff_score < low)

points(outliers.up$mean_score, outliers.up$diff_score, col = "blue")
points(outliers.lo$mean_score, outliers.lo$diff_score, col = "chocolate3")

dat$PID <- factor(dat$PID)

print(rmcorr(PID,mean_score, diff_score, dataset = dat))

}

else if(plot_option == FALSE) {

param <- c(low, mean, upper)

return(param)
}

else(NULL)
}

