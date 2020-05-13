demo_elastic <- function(training, testing) {
  
  train <- training %>% mutate(Cohort = "Train")
  test <- testing %>% mutate(Cohort = "Validation")
  
  dat <- bind_rows(train, test)
  
  dat$Cohort <- factor(dat$Cohort, levels = c("Train", "Validation"))
  
  age <- ggplot(dat, aes(x = Cohort, y = Age, fill = Cohort)) + geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "gainsboro") + theme_bw() + theme(legend.position = "none") +
    border("black") + theme(axis.ticks.length=unit(-0.25, "cm"),
                            axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                            axis.title.x = element_text(size = 14, face = "bold"),
                            axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label.x = 1.5)
  
  app <- ggplot(dat, aes(x = Cohort, y = App, fill = Cohort)) + geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "gainsboro") + theme_bw() + theme(legend.position = "none") +
    border("black") + theme(axis.ticks.length=unit(-0.25, "cm"),
                            axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                            axis.title.x = element_text(size = 14, face = "bold"),
                            axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label.x = 1.5)
  
  written <- ggplot(dat, aes(x = Cohort, y = Written, fill = Cohort)) + geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "gainsboro") + theme_bw() + theme(legend.position = "none") +
    border("black") + theme(axis.ticks.length=unit(-0.25, "cm"),
                            axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                            axis.title.x = element_text(size = 14, face = "bold"),
                            axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
    stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), label.x = 1.5)
  
  
  ggarrange(age, app, written, ncol = 3, nrow = 1)
  
  
}
