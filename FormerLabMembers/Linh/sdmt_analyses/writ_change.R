#-------------------------------------------------------
#Script to analyze clinically significant threshold
# in written SDMT 
#Linh Pham
#3/4/20
#--------------------------------------------------------

# set working directory: set this to be whatever location the 'sdmt' folder is located in
# for example, if you place this on your desktop, then set it to "...My Desktop/sdmt"
# setwd("h:/My Documents/sdmt_analyses")

# Load necessary packages
library(dplyr)
library(lubridate)
library(readr)
library(ggpubr)
library(ggplot2)
library(wesanderson)
library(tidyr)

# filter for trials that are greater than 0 days apart but less than 190 days apart
sixm <- read_csv("final_data/fig10a.csv")

# for each person, keep the first time that they had that days difference which fell into 
# the criteria of the sixm data. then calssify if the score difference is less than or greater than 4 points
score_change <- read_csv("final_data/fig10b.csv")

# calculate median in those below and above threshold and print
# the wilcoxon comparison values
median(subset(score_change, thresh == "<= -4")$diff_combi)
median(subset(score_change, thresh == "> -4")$diff_combi)
wilcox.test(score_change$diff_combi~score_change$thresh)

median(subset(score_change, thresh == "<= -4")$diff_neurex)
median(subset(score_change, thresh == "> -4")$diff_neurex)
wilcox.test(score_change$diff_neurex~score_change$thresh)

# calculate the new threshold based on IQR instead of the 4 points threshold that's currently being used
th <- quantile(score_change$diff_score)[2] - 1.5*IQR(score_change$diff_score)

# print out that threshold, which should be 13
print(th)

# new data frame that reclassifies if people fall below 13 or above 13 in their written score changes
sd2_change <- sixm %>% distinct(PID, .keep_all = TRUE) %>% 
              mutate(thresh = ifelse(diff_score <= -13, "<= -13", ">-13"))
              
# calculate median in those below and above threshold and print
# the wilcoxon comparison values                
median(subset(sd2_change, thresh == "<= -13")$diff_combi)
median(subset(sd2_change, thresh == ">-13")$diff_combi)
wilcox.test(sd2_change$diff_combi~sd2_change$thresh)

median(subset(sd2_change, thresh == "<= -13")$diff_neurex)
median(subset(sd2_change, thresh == ">-13")$diff_neurex)
wilcox.test(score_change$diff_neurex~score_change$thresh)

# plotting histogram and boxplots of the score differences in 6 months. this is 
# figure 10A in the paper
png(filename = "sdmt_figures/score_change.png", type = "cairo",
    units = "in", width = 6, height = 4, pointsize = 10, res = 400)

layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(score_change$diff_score , horizontal=TRUE , ylim=c(-22,11), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F)
par(mar=c(4, 3.1, 1.1, 2.1))
hist(score_change$diff_score , breaks=40 , col=rgb(0.2,0.8,0.5,0.5) , border=F , main="" , xlab="Score Differences in 6 Months", 
     xlim=c(-22,11), tck = 0.02)
abline(v = median(score_change$diff_score), col = "black", lwd = 2)
abline(v = -4, col = "black", lty = "dotted", lwd = 2)
abline(v = th, col = "brown", lty = "dashed", lwd = 2)

dev.off()

# plotting disability score differences based on the threshold of 4 points.
# this is figure 10B upper panel in the paper.
png(filename = "sdmt_figures/scale_diff.png", type = "cairo",
    units = "in", width = 9, height = 6, pointsize = 10, res = 400)

ed_4 <- ggplot(score_change, aes(x = thresh, y = diff_edss, fill = thresh)) + geom_violin(trim = FALSE) +
geom_boxplot(width = 0.1, fill = "white") +
scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest1")) +
theme_bw() + theme(legend.position = "none", axis.ticks.length=unit(-0.25, "cm"),
                   axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                   axis.title.x = element_text(size = 14, face = "bold"),
                   axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) 

combi_4 <- ggplot(score_change, aes(x = thresh, y = diff_combi, fill = thresh)) +
  geom_violin(trim = FALSE) + geom_boxplot(width=0.1, fill = "white") + 
  scale_fill_manual(values=wes_palette(n=2, name="Royal1")) +
  theme_bw() + theme(legend.position = "none", axis.ticks.length=unit(-0.25, "cm"),
                     axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                     axis.title.x = element_text(size = 14, face = "bold"),
                     axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) 

neur_4 <- ggplot(score_change, aes(x = thresh, y = diff_neurex, fill = thresh)) +
  geom_violin(trim = FALSE) + geom_boxplot(width=0.1, fill = "white") + 
  scale_fill_manual(values=wes_palette(n=2, name="Cavalcanti1")) +
  theme_bw() + theme(legend.position = "none", axis.ticks.length=unit(-0.25, "cm"),
                     axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                     axis.title.x = element_text(size = 14, face = "bold"),
                     axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))

# plotting disability score differences based on the threshold of 13 points.
# this is figure 10B lower panel in the paper.
ed_13 <- ggplot(sd2_change, aes(x = thresh, y = diff_edss, fill = thresh)) + geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest1")) +
  theme_bw() + theme(legend.position = "none", axis.ticks.length=unit(-0.25, "cm"),
                     axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                     axis.title.x = element_text(size = 14, face = "bold"),
                     axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) 

combi_13 <- ggplot(sd2_change, aes(x = thresh, y = diff_combi, fill = thresh)) +
  geom_violin(trim = FALSE) + geom_boxplot(width=0.1, fill = "white") + 
  scale_fill_manual(values=wes_palette(n=2, name="Royal1")) +
  theme_bw() + theme(legend.position = "none", axis.ticks.length=unit(-0.25, "cm"),
                     axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                     axis.title.x = element_text(size = 14, face = "bold"),
                     axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) 

neur_13 <- ggplot(sd2_change, aes(x = thresh, y = diff_neurex, fill = thresh)) +
  geom_violin(trim = FALSE) + geom_boxplot(width=0.1, fill = "white") + 
  scale_fill_manual(values=wes_palette(n=2, name="Cavalcanti1")) +
  theme_bw() + theme(legend.position = "none", axis.ticks.length=unit(-0.25, "cm"),
                     axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                     axis.title.x = element_text(size = 14, face = "bold"),
                     axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) 

ggarrange(ed_4, combi_4, neur_4, ed_13, combi_13, neur_13, ncol = 3, nrow = 2)

dev.off()
