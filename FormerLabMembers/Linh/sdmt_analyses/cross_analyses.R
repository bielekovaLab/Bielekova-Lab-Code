#-------------------------------------------------------
#Script to analyze symbol test cross-sectional data
#Linh Pham
#3/4/20
#--------------------------------------------------------

# set working directory: set this to be whatever location the 'sdmt' folder is located in
# for example, if you place this on your desktop, then set it to "...My Desktop/sdmt"
# setwd("/Users/linh/Desktop/sdmt_analyses-master")

# Load necessary packages
library(readr)
library(dplyr)
library(lubridate)  
library(ggplot2)
library(ggpubr)
library(car)
library(nlme)
library(caret)
library(glmnet)
library(Cairo)
library(rmcorr)
library(reshape)
library(DescTools)
library(ggcorrplot)

# Reading in all functions needed for cross-sectional analyses
pathnames <- list.files(pattern="[.]R$", 
             path="funcs_cross", 
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

#---------------------------------------------------------------------------------
#----------------------Begin Working with Cross-Sectional Data--------------------
#---------------------------------------------------------------------------------

# if you would like to see the resulting plots simply within R,
# run the line that would generate the plots and clear it by dev.off()
# but don't run the lines that start with png...
# to see the specifics of how each function works, simply go into the
# funcs_cross folder and find the function with the same name

#----------------------------------------------------------------------------------
#--------Obtaining the minimum time needed for reliable app SDMT scores------------
#----------------------------------------------------------------------------------

# get raw cross-sectional data that still maintains all guessess in each trial
# this dataset will be used in determining the minimum amount of time needed
# to obtain reliable app SDMT scores

split.cross.raw <- datcross(dat.combo)

# get a printout for the number of HV and MS in split.cross.raw
paste("split.cross.raw has", length(unique(subset(split.cross.raw, Diagnosis == "HV")$PID)), "HV and",
      length(unique(subset(split.cross.raw, Diagnosis == "MS")$PID)),"MS")

# Correlation of trial 1 v trial 2 over the test progression
# this output corresponds to figure 2B-C in the paper
png(filename = "results_cross/allCorPlots.png", type = "cairo",
    units = "in", width = 8, height = 4, pointsize = 10, res = 400)
allCorPlots(split.cross.raw)
dev.off()

#  Comparing 75sec and 90sec speed in MS and HV. This corresponds to Figure 2C in the paper. 
png(filename = "results_cross/compsec.png", type = "cairo",
    units = "in", width = 5, height = 5, pointsize = 10, res = 400)
compsec(split.cross.raw)
dev.off()

#----------------------------------------------------------------------------------
#-------- Outlier Adjustment and Cross-Sectional Data Generation ------------------
#----------------------------------------------------------------------------------

# Bland-altman plot of of ALL test sittings that have 2 app trials.
# This is Figure 3A in the paper
png(filename = "results_cross/bamix_orig.png", type = "cairo",
    units = "in", width = 4, height = 4, pointsize = 10, res = 400)
mix_ba1v2(dat.combo,plot_option = TRUE)
dev.off()

# print Bland-altman limits of agreement and mean difference
mix_ba1v2(dat.combo, plot_option = FALSE)

# Bland-altman with outliers adjusted dataset. Limits of agreement and mean lines
# came from the lines set by the raw data set (unadjusted)
# This is figure 3C in the paper
png(filename = "results_cross/bamix_adj.png", type = "cairo",
    units = "in", width = 4, height = 4, pointsize = 10, res = 400)
outlier.adj(dat.combo, option = TRUE)
dev.off()

# Obtain full data set with outliers adjusted 
# see methods section of the paper or function outlier_adj for more details
full.orm <- outlier.adj(dat.combo, option = FALSE)

# create working cross-sectional dataset with outliers adjusted
# this is the cross-sectional dataset that will be used for all 
# remaining cross-sectional analyses
cross <- full.orm %>% ungroup() %>% group_by(PID) %>% filter(Date == min(Date)) %>%
  group_by(PID, Date) %>% filter(Time.Vals == min(Time.Vals))

# get a printout for the number of HV and MS in cross-sectional data
paste("cross has", length(unique(subset(cross, Diagnosis == "HV")$PID)), "HV and",
      length(unique(subset(cross, Diagnosis == "MS")$PID)),"MS")

#----------------------------------------------------------------------------------
#----------------------Validity and Clinical Association Analyses------------------
#----------------------------------------------------------------------------------

# Comparing app scores between MS and HV. This is Figure 4A in the paper. 
png(filename = "results_cross/hvms.png", type = "cairo",
    units = "in", width = 3, height = 3, pointsize = 10, res = 400)
hvmsbpx(cross)
dev.off()

# Obtain the concordance plot between app and written scores. 
# This is Figure 4B in the paper.

# grabbing data needed for fig4B-C analyses
fig4B <- read_csv("final_data/fig4B.csv")
fig4C <- read_csv("final_data/fig4C.csv")

png(filename = "results_cross/corAvW_orig.png", type = "cairo",
     units = "in", width = 4, height = 4, pointsize = 10, res = 400)
corAvW(cross, "orig")
dev.off()

# Obtain the Bland-Altman plot between the written and app scores.
# This is Figure 4C in the paper. 

png(filename = "results_cross/baAvW.png", type = "cairo",
    units = "in", width = 4, height = 4, pointsize = 10, res = 400)
baAvW(cross,option = TRUE)
dev.off()

# Obtain the concordance plot between app and written scores,
# with the app scores now increased by 7.8 points to 
# account for the difference betwen written and app scores in general
# Figure 4D in the paper
png(filename = "results_cross/corAvW_adj.png", type = "cairo",
    units = "in", width = 4, height = 4, pointsize = 10, res = 400)
corAvW(cross, "adj")
dev.off()

# Obtain correlation matrix between the PASAT-3, written, and app scores
# with clinical measurements. this is figure 5 in the paper. 
fig5 <- read_csv("final_data/fig5.csv")

png(filename = "results_cross/matrix.png", type = "cairo",
    units = "in", width = 10, height = 4, pointsize = 10, res = 400)
make_matrix(cross, plot_option = TRUE)
dev.off()

# print out how many HV and MS individuals are included in the
# correlation matrix analyses. save the merged data(smartphone + all clinical + mri)
# to an object as well
cross_clin <- make_matrix(cross, plot_option = FALSE)

# run the following line to make sure seed is set consistently
# across all R versions
RNGkind(sample.kind = "Rounding")

# elastic net regression analyses start below this line
# first step is to make the training dataset from
# the individuals who were used in the correlation matrix analyses
trainset <- read_csv("final_data/trainset.csv")

# get a printout for the number of HV and MS in the training set
paste("trainset has", nrow(subset(trainset, Diagnosis == "HV")), "HV and",
      nrow(subset(trainset, Diagnosis == "MS")),"MS")

# create a validation dataset from the individuals 
# who were used in the correlation matrix analyses
testset <- read_csv("final_data/testset.csv")

# get a printout for the number of HV and MS in the training set
paste("testset has", nrow(subset(testset, Diagnosis == "HV")), "HV and",
      nrow(subset(testset, Diagnosis == "MS")),"MS")

# produce elastic net results for the written SDMT. figure 6A in the paper
png(filename = "results_cross/written_elastic.png", type = "cairo",
    units = "in", width = 4, height = 1.8, pointsize = 10, res = 400)
elastic_results(trainset,testset,"Written", "neurex")
dev.off()

# produce elastic net results for app SDMT. figure 6B in the paper
png(filename = "results_cross/app_elastic.png", type = "cairo",
    units = "in", width = 4, height = 1.8, pointsize = 10, res = 400)
elastic_results(trainset,testset,"App","neurex")
dev.off()

# produce app smdt elastic net results with tapping as surrogate for
# hand mobility functions. figure 7A in the paper
png(filename = "results_cross/app_tap_elastic.png", type = "cairo",
    units = "in", width = 4, height = 1.8, pointsize = 10, res = 400)
elastic_results(trainset,testset,"App", "phone")
dev.off()

# getting sdmt-tap elastic net result coefficients that are UNSTANDARDIZED (output not scaled)
elastic_unscaled(trainset,testset, "App", "phone")

# produce app elastic net results without tapping as surrogate for 
# hand mobility functions. figure 7B in the paper. 
png(filename = "results_cross/app_no_tap_elastic.png", type = "cairo",
    units = "in", width = 4, height = 1.8, pointsize = 10, res = 400)
elastic_results(trainset,testset,"App", "no tap")
dev.off()

#-------------------------------------------------------------------#
#------------------------Make Supplementary Figures-----------------#
#-------------------------------------------------------------------#

# suppl.2-4: correlation plots of all the variables that are in the correlation matrix
png(filename = "suppl_figs/pasat_corr.png", type = "cairo",
    units = "in", width = 10, height = 10, pointsize = 10, res = 400)
plot_matrix(cross, "PASAT")
dev.off()

png(filename = "suppl_figs/app_corr.png", type = "cairo",
    units = "in", width = 10, height = 10, pointsize = 10, res = 400)
plot_matrix(cross, "App")
dev.off()

png(filename = "suppl_figs/writ_corr.png", type = "cairo",
    units = "in", width = 10, height = 10, pointsize = 10, res = 400)
plot_matrix(cross, "Written")
dev.off()


# suppl.6: distribution of age and app/written sdmt scores in the
# training and validation cohorts

png(filename = "suppl_figs/demo_elastic.png", type = "cairo",
    units = "in", width = 10, height = 3, pointsize = 10, res = 400)
demo_elastic(trainset,testset)
dev.off()

# suppl. 7A. Diagnostic plots for written SDMT elastic net
png(filename = "suppl_figs/diag_writ_neurex.png", type = "cairo",
    units = "in", width = 4.5, height = 8, pointsize = 10, res = 400)
elastic_diag(trainset,testset,"Written", "neurex")
dev.off()

# suppl. 7B. Diagnostic plots for app SDMT elastic net
png(filename = "suppl_figs/diag_app_neurex.png", type = "cairo",
    units = "in", width = 4.5, height = 8, pointsize = 10, res = 400)
elastic_diag(trainset,testset,"App","neurex")
dev.off()

# suppl. 8A. Diagnostic plots for app SDMT with tapping elastic net
png(filename = "suppl_figs/diag_app_tap.png", type = "cairo",
    units = "in", width = 4.5, height = 8, pointsize = 10, res = 400)
elastic_diag(trainset,testset,"App", "phone")
dev.off()

# suppl. 8B. Diagnostic plots for app SDMT without tapping elastic net
png(filename = "suppl_figs/diag_app_notap.png", type = "cairo",
    units = "in", width = 4.5, height = 8, pointsize = 10, res = 400)
elastic_diag(trainset,testset,"App", "no tap")
dev.off()

# suppl. 12. Elastic net results with upper motor strength included
elastic_results_all(trainset,testset,"App","neurex")

# suppl. 13. Written SDMT comparison between HV and MS
png(filename = "suppl_figs/bpWvW.png", type = "cairo",
    units = "in", width = 4, height = 4, pointsize = 10, res = 400)
bpWvW(cross)
dev.off()
