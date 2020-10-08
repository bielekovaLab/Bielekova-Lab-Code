test <- full.orm %>% ungroup() %>% group_by(PID, Date) %>% summarise(score = mean(score_p_sec)*90) %>%
        ungroup() %>% group_by(PID) %>% mutate(sitting = rank(Date)) %>%
        filter(sitting <= 2) %>% filter(n() == 2) %>% 
        mutate(days_diff = Date - min(Date)) %>% 
        mutate(rm = ifelse(max(days_diff) < 160, "T", "F")) %>%
        filter(rm == "F")

t1 <- test %>% filter(sitting == 1)
t2 <- test %>% filter(sitting == 2)

plot(t2$score~t1$score, xlab = "Month 0 Score", 
     ylab = "~Month 6 Score", main = "App SDMT Scores")

cor.test(t1$score, t2$score, method = "spearman")
summary(lm(t2$score~t1$score))

binding <- data.frame(t1$score, t2$score)

print(psych::corr.test(binding, use = "pairwise", 
                 method = "spearman", adjust = "none"), short = FALSE)

test <- test %>% select(-rm)
write_csv(test, "extra_analyses/request_0602.csv")


writ_test <- read_csv("final_data/fig10a.csv") %>% 
  group_by(PID) %>% filter(diff > 160) %>% distinct(PID, .keep_all = TRUE) %>%
  mutate(t1 = Paper_Correct + diff_score)
  

binding <- data.frame(writ_test$t1, writ_test$Paper_Correct)

print(psych::corr.test(binding, use = "pairwise", 
                       method = "spearman", adjust = "none"), short = FALSE)
