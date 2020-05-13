#-------------------------------------------------------
#Script to analyze symbol test longitudinal data
#Linh Pham
#3/4/20
#--------------------------------------------------------

# set working directory: set this to be whatever location the 'sdmt' folder is located in
# for example, if you place this on your desktop, then set it to "...My Desktop/sdmt"
setwd("h:/My Documents/sdmt_analyses")

# Load necessary packages
library(readr)
library(dplyr)
library(lubridate)  
library(ggplot2)
library(ggpubr)
library(car)
library(nlme)
library(Cairo)
library(tidyr)
library(sjPlot)
library(lme4)
library(lattice)
library(zoo)
library(purrr)
library(gridExtra)


# Reading in all functions needed for longitudinal analyses
pathnames <- list.files(pattern="[.]R$", 
                        path="funcs_long", 
                        full.names=TRUE);
                        sapply(pathnames, FUN=source)

# Reading in data:

# sec75v1 is app sdmt data from 75seconds trial, app version 1.
# sec90pre1 is app sdmt data from 90seconds trial, before app version 1.
# sec90_1 is app sdmt data from 90seconds trial, starting from app version 1. 
sec75v1 <- read.csv("data_clean/75secV1.csv")
sec90pre1 <- read.csv("data_clean/90secPreV1.csv")
sec90_1 <- read.csv("data_clean/90secV1.csv")

# ------------------Begin Data Cleaning----------------------------------------
# 1. Combine the 90sec data from 2 different app versions together
# 2. Create a new column called 'score_p_sec': value here is #correct/total 
# time spent on test (in seconds)
# 3. Call this new combined data 'sec90'. Do the same procedure to find 
# correct/total time spent on test per trial for the 75 seconds data. 
# Difficulty for 75sec is filtered to be 1 because that's the full contrast trial
# (we've developed a version of the app SDMT with medium gray/red contrasted symbols,
# and the difficulty for this version is either 2 or 3. we're not analyzing
# these difficulty levels for now). 
#---------------------------------------------------------------------------------

sec90 <- rbind(sec90pre1,sec90_1) %>% group_by(PID,Date,Time) %>% 
         mutate(score_p_sec = Correct/90) %>% ungroup() 

sec75 <- sec75v1 %>% group_by(PID, Date, Time) %>% 
         mutate(score_p_sec = Correct/75) %>% 
         filter(Difficulty == 1) %>% ungroup()

# combining the two datasets together, with a column called "Groups'
# which specify the origin of each row.
dat.combo <- bind_rows("90sec" = sec90,"75sec" = sec75, .id = "Groups")

# making sure the date format is correct for working in R
dat.combo$Date <- ymd(dat.combo$Date)

#----------------------------------------------------------------------------------
#-------- Outlier Adjustment and Longitudinal Data Generation ---------------------
#----------------------------------------------------------------------------------

# Obtain full data set with outliers adjusted 
# see methods section of the paper or function outlier_adj for more details
full.orm <- outlier.adj(dat.combo, option = FALSE)

# Obtaining longitudinal data (individuals with >= 20 trials)
all.long.rm <- all_long(dat.combo, "no age")

# get the same longitudinal data set but with age included for demographics purposes
all.long.demo <- all_long(dat.combo, "with age")

#---------------------------------------------------------------------
#-------------------Begin Longitudinal Analyses-----------------------
#---------------------------------------------------------------------

#------------------------------------------------------
# Finding learning stop points in longitudinal data
#------------------------------------------------------

# Defining the non-linear regression function for identifying 
# the change point, c

fx <- function(x,c,b0,b1){
  out <- b0 + as.numeric(x < c)*b1*x + as.numeric(x>=c)*b1*c
  return(out)
}

# Changing column headers, average daily scores, and creating a column x which
# signifies sitting number in  the longitudinal data so that it will work with the 
# non-linear regression function 
transformed <- long.transform(all.long.rm)

# Plot the longitudinal data and indicate change point 
# based on the non-linear regression function
# this is figure 8 in the paper
png(filename = "results_long/all_cpt.png", type = "cairo",
    units = "in", width = 9, height = 6, pointsize = 10, res = 400)
plot_cpt(print_out = TRUE)
dev.off()

# get a printout of all individuals spearman correlation after the learning period
# and the p-values associated with those correlations
plot_cpt(print_out = FALSE)

# ----------------------------------------- #
# ICC Analyses for Longitudinal Patients    #
# ----------------------------------------- #

# obtain all of the longitudinal data (with scores averaged), including the 
# period that the patients are still learning
pre_learn <- all.long.rm %>% group_by(PID, Date, Time) %>%
  ungroup() %>% group_by(PID, Date) %>%
  summarise(score = mean(Correct)*90) %>%
  ungroup() %>% group_by(PID) %>% mutate(sitting = rank(Date))

# obtain the longitudinal data (with scores averaged), NOT including the 
# period that the patients are still learning (ONLY post-learning is included)
post_learn <- all.long.rm %>% group_by(PID, Date, Time) %>%
        ungroup() %>% group_by(PID, Date) %>%
        summarise(score = mean(Correct)*90) %>%
        ungroup() %>% group_by(PID) %>% mutate(sitting = rank(Date)) %>%
  filter(PID == "AAO02" & sitting > 6 |
      PID == "ABX699" & sitting > 3 |
      PID == "ABX711" & sitting > 10 |
      PID == "ACW072" & sitting > 4  |
      PID == "ACW204" & sitting > 12 |
      PID == "ACW254" & sitting > 6 |
      PID == "ACW269" & sitting > 3 |
      PID == "ACW282" & sitting > 5  |
      PID == "ACW365" & sitting > 5  |
      PID == "ACW406" & sitting > 13|
      PID == "ACW407" & sitting > 8 |
      PID == "ACW484" & sitting > 4 |
      PID == "ACW570" & sitting > 12|
      PID == "ABX710" & sitting > 15)

# mixed model to get ICC values with learning included and
# without learning included
m1 <- lmer(score ~ sitting + (1|PID), data = pre_learn)
m2 <- lmer(score ~ sitting + (1|PID), data = post_learn)

# ICC values of the two different mixed models
tab_model(m1)
tab_model(m2)

# begin plotting ICC analyses results
# here, we wanted to plot both the app SDMT data AND app SDMT scores
# converted into written scores, using the mean difference
# between the two methods, obtained in the cross-sectional 
# analysis (see figure 4A in the paper)
thresh.mean <- 7.777717

# preparing the data for plotting
# create a new dataframe called pre, which is simply
# all of the longitudinal data in the pre_learn dataframe + a new column that 
# indicates this comes from the app data
pre <- pre_learn %>% mutate(type = "app")

# create a new dataframe called pre.adj, which is
# all of the longitudinal data in the pre_learn dataframe with
# the conversion to written score using the mean difference (thresh.mean) and a 
# new column that indicates this column is the written scores converted from app scores
pre.adj <- pre_learn %>% mutate(score = score + thresh.mean,
                                type = "converted")

# combine these two previous dataframes together
pre_learn <- bind_rows(pre, pre.adj)

# do the same adjustments for the post_learning data
post <- post_learn %>% mutate(type = "app")
post.adj <- post_learn %>% mutate(score = score + thresh.mean,
                                  type = "converted")

# because the original post_learning data didn't include ACW245 and ACW559
# we're adding them here but keeping their data as NA. on the final plot,
# this will show them as blank data
post_learn_plot <- bind_rows(post, post.adj) %>% ungroup() %>%
                   add_row(PID = c("ACW245", "ACW245", "ACW559", "ACW559"),
                           Date = c(NA,NA,NA,NA),
                           score = c(NA, NA,NA,NA),
                           sitting = c(NA, NA,NA,NA),
                           type = c("app", "converted", "app", "converted"))


# plot the longitudinal data as a strip chart and indicate the ICC values
# associated with all of the longitudinal data or with just the post-learning
# data. this is figure 9 in the paper. 
png(filename = "results_long/icc.png", type = "cairo",
    units = "in", width = 8, height = 6, pointsize = 10, res = 400)

p <- ggplot(pre_learn, aes(x = PID, y = score, color = type, shape = type)) + 
  scale_shape_manual(values=c(1,2))+ 
geom_jitter(position=position_dodge(0.7)) + theme_minimal() + 
  theme(legend.position = "none")

q <- ggplot(post_learn_plot, aes(x = PID, y = score, color = type, shape = type)) + 
  scale_shape_manual(values=c(1,2))+
  geom_jitter(position=position_dodge(0.6)) + theme_minimal() + theme(legend.position = "none")

ggarrange(p, q, nrow = 2)

dev.off()

#------------------------------------------------------#
#--------Finding clinically Significant Threshold------#
#-------------------in app SDMT------------------------#
#------------------------------------------------------#

# find unique patient ids in the post-learning data
ids <- unique(post_learn$PID)

# find the means score of 1 sitting, 2 sittings,
# 3 sittings, 4 sittings, for all longitudinal, post-learning patients
all_means <- bind_rows(lapply(1:length(ids),function(x) group_means(post_learn, ids[x])))

# plotting a demonstration of how the grouping average is done. this is figure 11A. 
png(filename = "results_long/means_diff_demo.png", type = "cairo",
    units = "in", width = 8, height = 2, pointsize = 10, res = 400)
demo_pt()
dev.off()

# find the differences of the mean scores for the grouped sittings
# eliminate the difference from the first sitting since it will be calculated as 0 
# by R
means_diff <- all_means %>% group_by(grouping, pid) %>% arrange(x) %>%
              mutate(diff = y - lag(y, default = first(y))) %>%
              filter(x != 1)

# plot the differences calculated in means_diff. this is figure 11B in the paper. 
png(filename = "results_long/means_diff_plot.png", type = "cairo",
    units = "in", width = 5, height = 4, pointsize = 10, res = 400)
plot_means()
dev.off()

#-----------------------------------------------------#
#------------Making Supplementary Figures-------------#
#------Associated with Longitudinal Data--------------#
#-----------------------------------------------------#

# Diagnostic plots for mixed bland-altman equations
# upper plots are for the equations that make the limits of agreement
# lower plots are for the equations that make the mean difference
# this is supplementary figure 1
png(filename = "suppl_figs/mix_ba_DiagPlots.png", type = "cairo",
    units = "in", width = 8, height = 8, pointsize = 10, res = 400)
mix_ba_diag(dat.combo)
dev.off()


# diagnostic plots of model. 
# residuals are clumping due to small amount of patients being analyzed
# this is supplementary figure 9
a <- plot(m1)
b <- plot(m2)
c <- qqmath(m1)
d <- qqmath(m2)

png(filename = "suppl_figs/mix_ICC_DiagPlots.png", type = "cairo",
    units = "in", width = 8, height = 8, pointsize = 10, res = 400)
grid.arrange(a, c, b, d, ncol = 2, nrow = 2)
dev.off()

# plot of longitudinal data based on days from first trial instead of sitting number
# this is supplementary figure 10
png(filename = "suppl_figs/long_by_days.png", type = "cairo",
    units = "in", width = 9, height = 6, pointsize = 10, res = 400)
delta.days(transformed)
dev.off()

# ICC analyses based on time diference instead of sitting number
# this is part of the legend in supplementary figure 11
pre_diff <- pre_learn %>% group_by(PID) %>% mutate(diff = Date - min(Date)) %>% filter(type == "app")
post_diff <- post_learn %>% group_by(PID) %>% mutate(diff = Date - min(Date))

m1_diff <- lmer(score ~ diff + (1|PID), data = pre_diff)
m2_diff <- lmer(score ~ diff + (1|PID), data = post_diff)

a <- plot(m1_diff)
b <- plot(m2_diff)
c <- qqmath(m1_diff)
d <- qqmath(m2_diff)

# diagnostic plots of ICC difference based on the days difference
# this is supplementary figure 11
png(filename = "suppl_figs/mix_diffICC_DiagPlots.png", type = "cairo",
    units = "in", width = 8, height = 8, pointsize = 10, res = 400)
grid.arrange(a, c, b, d, ncol = 2, nrow = 2)
dev.off()

# ICC results
tab_model(m1_diff) # ICC = 0.88
tab_model(m2_diff) # ICC = 0.91
