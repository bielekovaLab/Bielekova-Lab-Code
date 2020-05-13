hvmsbpx <- function(data) { 
app_mean <- data %>% group_by(PID, Date, Diagnosis, Gender, Age) %>% 
            summarise(score_p_sec = mean(score_p_sec)*90)

yrange <- range(app_mean$score_p_sec)

label.pos <- max(yrange)

diag.statement <- paste("HV median is", median(subset(app_mean, Diagnosis == "HV")$score_p_sec),
                        ",HV range is", range(subset(app_mean, Diagnosis == "HV")$score_p_sec)[1], "to", 
                        range(subset(app_mean, Diagnosis == "HV")$score_p_sec)[2],
                        ",MS median is", median(subset(app_mean, Diagnosis == "MS")$score_p_sec),
                        "and MS range is", range(subset(app_mean, Diagnosis == "MS")$score_p_sec)[1], "to",
                        range(subset(app_mean, Diagnosis == "MS")$score_p_sec)[2])

print(diag.statement)
print(paste("HV mean is", mean(subset(app_mean, Diagnosis == "HV")$score_p_sec)))
print(paste("MS mean is", mean(subset(app_mean, Diagnosis == "MS")$score_p_sec)))


p <- ggplot(app_mean, aes(x = Diagnosis, y = score_p_sec, fill = Diagnosis)) + geom_violin(trim = FALSE) + 
     scale_fill_manual(values = c("dodgerblue2", "chocolate3")) +
     geom_boxplot(width = 0.1, fill = "gainsboro") + theme_bw() + theme(legend.position = "none") +
     border("black") + theme(axis.ticks.length=unit(-0.25, "cm"),
                          axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                          axis.title.x = element_text(size = 14, face = "bold"),
                          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) +
     ylim(-10, 150)

print(p)

print(paste("number of HV:", length(unique(subset(app_mean, Diagnosis == "HV")$PID))))
print(paste("number of MS:", length(unique(subset(app_mean, Diagnosis == "MS")$PID))))

print(wilcox.test(score_p_sec~Diagnosis, data = app_mean))
} 
