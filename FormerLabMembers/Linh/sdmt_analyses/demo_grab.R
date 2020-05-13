#----------------------------------
# script to grab demographics info
# 03/05/2020
#---------------------------------

# set working directory
# setwd("h:/My Documents/sdmt_analyses")

# grab the functions that are needed for running the edss/neurex matching process
source("funcs_cross/match_clin.R")
source("funcs_cross/app_match_mean.R")

# load necessary packages
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)

# read in data to obtain demographics measurement

# all cross-sectional smartphone data
cross_all <- read_csv("demographics/cross_all.csv") %>% dplyr::select(PID, Age, Gender, Date) %>%
             distinct(PID, .keep_all = TRUE)

# cross-sectional smartphone data with ALL clinical values
cross_clin <- read_csv("demographics/cross_clin.csv") %>% dplyr::select(PID, Age, Gender, Date) %>%
              distinct(PID, .keep_all = TRUE)

# longitudinal smartphone data
long <- read_csv("demographics/long.csv") %>% group_by(PID) %>% filter(Date == min(Date)) %>%
        ungroup() %>% dplyr::select(PID, Date, Age) %>% distinct(PID, .keep_all = TRUE)


# reading in written threshold demographics data
writ_tbl <- read_csv("demographics/written_threshold_demographics.csv")

# reading neurex and edss data
clinical <- read_csv("data_raw/neurex.raw.csv") %>% dplyr::select(patientcode, date, edss_score_calculated,
                                                                  neurex_total) %>%
            dplyr::rename("PID" = "patientcode", 
                          "Date" = "date",
                          "edss" = "edss_score_calculated",
                          "neurex" = "neurex_total")

# reformat date in the clinical data
clinical$Date <- mdy(clinical$Date)

# reading in gender data
gender <- read_csv("data_clean/handedness.csv") %>% dplyr::select(PatientCode, Gender) %>% dplyr::rename("PID" = "PatientCode")

# reading in detailed diagnosis data
diag <- read_csv("data_clean/diag_detailed.csv")

#-----------------------------------------------------------
#-------------Demographics for cross-sectional ALL----------
#-----------------------------------------------------------
crossAll_tbl <- read_csv("demographics/crossAll_tbl.csv")
crossAll_tbl$diagnosis <- factor(crossAll_tbl$diagnosis, levels = c("PP-MS", "SP-MS", "RR-MS", "HV"))

# print age information
crossAll_tbl %>% group_by(diagnosis) %>% summarise(age = mean(Age),
                                                   sd = sd(Age),
                                                   min = min(Age),
                                                   max = max(Age))
# print gender information
crossAll_tbl %>% group_by(diagnosis, Gender) %>% tally() %>% ungroup() %>% group_by(diagnosis) %>% 
mutate(total = sum(n)) %>% mutate(percent = round((n/total)*100))


#----------------------------------------------------------
#------Demographics for cross-sectional with clinical------
#----------------------------------------------------------
crossClin_tbl <- read_csv("demographics/crossClin_tbl.csv")
crossClin_tbl$diagnosis <- factor(crossClin_tbl$diagnosis, levels = c("PP-MS", "SP-MS", "RR-MS", "HV"))

# print age information
crossClin_tbl %>% group_by(diagnosis) %>% summarise(age = mean(Age),
                                                   sd = sd(Age),
                                                   min = min(Age),
                                                   max = max(Age))
# print gender information
crossClin_tbl %>% group_by(diagnosis, Gender) %>% tally() %>% ungroup() %>% group_by(diagnosis) %>% 
  mutate(total = sum(n)) %>% mutate(percent = round((n/total)*100))

# print neurex info
crossClin_tbl %>% group_by(diagnosis)  %>% summarise(avg = mean(neurex),
                                                    sd = sd(neurex),
                                                    min = min(neurex),
                                                    max = max(neurex))
# print edss info
crossClin_tbl %>% group_by(diagnosis) %>% summarise(avg = mean(edss),
                                                    sd = sd(edss),
                                                    min = min(edss),
                                                    max = max(edss))

#-----------------------------------------------------------------
#-------Demographics for written threshold analysis---------------
#-----------------------------------------------------------------
writ_tbl$diagnosis <- factor(writ_tbl$diagnosis, levels = c("PP-MS", "SP-MS", "RR-MS", "HV"))

# print age information
writ_tbl %>% group_by(diagnosis) %>% summarise(age = mean(Age),
                                                    sd = sd(Age),
                                                    min = min(Age),
                                                    max = max(Age))
# print gender information
writ_tbl %>% group_by(diagnosis, Gender) %>% tally() %>% ungroup() %>% group_by(diagnosis) %>% 
  mutate(total = sum(n)) %>% mutate(percent = round((n/total)*100))

# print neurex info
writ_tbl %>% group_by(diagnosis) %>% summarise(avg = mean(neurex),
                                                    sd = sd(neurex),
                                                    min = min(neurex),
                                                    max = max(neurex))
# print edss info
writ_tbl %>% group_by(diagnosis) %>% summarise(avg = mean(edss),
                                                    sd = sd(edss),
                                                    min = min(edss),
                                                    max = max(edss))

#-------------------------------------------------------------
#--------Demographics for longitudinal analysis---------------
#-------------------------------------------------------------
long_tbl <- read_csv("demographics/long_tbl.csv")
long_tbl$diagnosis <- factor(long_tbl$diagnosis, levels = c("PP-MS", "SP-MS", "RR-MS", "HV"))

# print age information
long_tbl %>% group_by(diagnosis) %>% summarise(age = mean(Age),
                                               sd = sd(Age),
                                               min = min(Age),
                                               max = max(Age))

# print gender information
long_tbl %>% group_by(diagnosis, Gender) %>% tally() %>% ungroup() %>% group_by(diagnosis) %>% 
  mutate(total = sum(n)) %>% mutate(percent = round((n/total)*100))


# merge to get neurex and edss info + add 1 more row for ACW282 since 
# its neurex and edss data fell outside of the matching threshold (2-3 days). it's actually 5 days for this person. 
long_fin <- app_match_mean(long_tbl, clinical) %>% dplyr::select(PID, diagnosis, edss, neurex) %>%
            add_row(PID = "ACW282", diagnosis = "RR-MS", edss = 3.5, neurex = 58.95555)

# print neurex info
long_fin %>% group_by(diagnosis) %>% summarise(avg = mean(neurex),
                                               sd = sd(neurex),
                                               min = min(neurex),
                                               max = max(neurex))
# print edss info
long_fin %>% group_by(diagnosis) %>% summarise(avg = mean(edss),
                                               sd = sd(edss),
                                               min = min(edss),
                                               max = max(edss))
