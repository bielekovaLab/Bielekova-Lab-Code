# R packages utilized
library(tidyverse)
library(ranger)
library(lubridate)

# Function to split a character representing a ratio in separate pieces
split_rat <- function(rat){
  split_list <- unlist(strsplit(rat,"[/]"))
  v1 <- first(split_list)
  v2 <- last(split_list)
  return(c(v1,v2))
}

# Function to take adjust a list of somamer ratios (ratio list) using healthy donar data (hd_data)
# and append them to the end of a dataset (data)
ratio_adjust_append <- function(data,ratio_list,hd_data){
  p <- length(ratio_list)
  for(i in 1:p){
    mark <- split_rat(ratio_list[i])
    mod1 <- adjust_function(mark[1],data=hd_data)
    mod2 <- adjust_function(mark[2],data=hd_data)
    out <- (log(data[[mark[1]]]) - predict(mod1,newdata=data)) - (log(data[[mark[2]]]) - predict(mod2,newdata=data))
    data[,ratio_list[i]] <- out
  }
  return(data)
}

# Function to create a model predicting a somamer (y_char) using healthy donor data (hd_soma)
# if it is in a list of somamers to adjust (age_to_adjust, sex_to_adjust)
adjust_function <- function(y_char,data=hd_soma,alpha_val=0.05,...){
  if(nchar(y_char)==17){data <- ratio_append(data,y_char,expo = FALSE)}
  y<-log(data[[y_char]])
  if(any(str_detect(y_char,c(age_to_adjust,sex_to_adjust)))){
    
    if(any(str_detect(y_char,age_to_adjust)) & any(str_detect(y_char,sex_to_adjust))){
      mod <- lm(y~age+gender,data)
    }else{
      if(any(str_detect(y_char,age_to_adjust))){
        mod <- lm(y~age,data)
      }
      if(any(str_detect(y_char,sex_to_adjust))){
        mod <- lm(y~gender,data)
      }
    }
  }else{
    mod <- lm(y~1)
  }
  return(mod)
}

# Load MS data and perform slight cleaning
soma_untreated <- read_csv("./inputs/CLEANED_soma_untreated_20190508.csv")
soma_untreated <- soma_untreated %>% 
  mutate(gender = factor(gender,levels = c("M","F"))) %>% 
  mutate(age_mri = as.numeric(mdy(mri_date)-mdy(dob))/365)

# Create MRI severity measure (needs to be inverted later)
mri_mod <- lm(volume_calculated_brainparenchymalfr~age_mri,data=soma_untreated)
soma_untreated$mri_severity <- soma_untreated$volume_calculated_brainparenchymalfr - predict(mri_mod,newdata=soma_untreated)

# Load in HD data
hd_soma <- read_csv("./inputs/CLEANED_HD_20181109.csv")

# Load in lists of somamers with HD relationships with age and gender
age_to_adjust <- readLines("./inputs/age_to_adjust_20190626.txt")
sex_to_adjust <- readLines("./inputs/sex_to_adjust_20190626.txt")

# Load in the final list of ratios for the three pipelines
msdss_baseline_ratios <- readLines("./inputs/msdss_baseline_ratios.txt")
msdss_followup_ratios <- readLines("./inputs/msdss_followup_ratios.txt")
mri_severity_baseline_ratios <- readLines("./inputs/mri_severity_baseline_ratios.txt")

# Create a list of all the ratios from the three pipelines
adj_ratios <- c(msdss_baseline_ratios,msdss_followup_ratios,mri_severity_baseline_ratios)

# Add these ratios to the master list, then split into training and validation
soma_untreated <- soma_untreated %>% ratio_adjust_append(adj_ratios,hd_data = hd_soma)
soma_train <- soma_untreated %>% filter(model_cohort=="training")
soma_test <- soma_untreated %>% filter(model_cohort=="validation")

# Create model of MS-DSS at baseline, save to outputs
msdss_baseline_data <- soma_train %>% 
  select(msdss,which(names(soma_train) %in% msdss_baseline_ratios))
nvars <- floor(3*sqrt(dim(msdss_baseline_data)[2]))
msdss_baseline_mod <- ranger(data=as.data.frame(msdss_baseline_data), num.trees=40000,mtry=nvars,importance='impurity',
              seed=223,dependent.variable.name="msdss",write.forest = TRUE)
saveRDS(msdss_baseline_mod,"./outputs/msdss_baseline_model.rds")

# Create model of MS-DSS at followup, save to outputs
msdss_followup_data <- soma_train %>% 
  select(msdss_last,which(names(soma_train) %in% msdss_followup_ratios))
nvars <- floor(3*sqrt(dim(msdss_followup_data)[2]))
msdss_followup_mod <- ranger(data=as.data.frame(msdss_followup_data), num.trees=40000,mtry=nvars,importance='impurity',
                             seed=223,dependent.variable.name="msdss_last",write.forest = TRUE)
saveRDS(msdss_followup_mod,"./outputs/msdss_followup_model.rds")

# Create model of MRI severity at baseline, save to outputs
mri_severity_baseline_data <- soma_train %>% 
  select(mri_severity,which(names(soma_train) %in% mri_severity_baseline_ratios)) %>% 
  filter(!is.na(mri_severity))
nvars <- floor(3*sqrt(dim(mri_severity_baseline_data)[2]))
mri_severity_baseline_mod <- ranger(data=as.data.frame(mri_severity_baseline_data), num.trees=40000,mtry=nvars,importance='impurity',
                             seed=223,dependent.variable.name="mri_severity",write.forest = TRUE)
saveRDS(mri_severity_baseline_mod,"./outputs/mri_severity_baseline_model.rds")

# Add OOB predictions to the training dataset
soma_train$soma_msdss_adjust <- msdss_baseline_mod$predictions
soma_train$soma_msdss_last_adjust <- msdss_followup_mod$predictions

# Round-about way to add predictions for MRI severity, since it is missing in some patients
soma_train_mri <- filter(soma_train,!is.na(mri_severity))
soma_train_mri$soma_mri_severity <- (-1)*mri_severity_baseline_mod$predictions
soma_train_mri <- select(soma_train_mri,sampleid,soma_mri_severity)
soma_train <- left_join(soma_train,soma_train_mri,by=c("sampleid"))

# Add model predictions to the validation cohort
soma_test$soma_msdss_adjust <- predict(msdss_baseline_mod,soma_test)$predictions
soma_test$soma_msdss_last_adjust <- predict(msdss_followup_mod,soma_test)$predictions
soma_test$soma_mri_severity <- predict(mri_severity_baseline_mod,soma_test)$predictions

# Invert the MRI severity measures
soma_test$soma_mri_severity <- (-1)*soma_test$soma_mri_severity
soma_test$mri_severity <- (-1)*soma_test$mri_severity
soma_train$mri_severity <- (-1)*soma_train$mri_severity

# Test the performance in the validation cohort
with(soma_test,cor.test(msdss,soma_msdss_adjust,method="spearman"))
with(soma_test,cor.test(msdss_last,soma_msdss_last_adjust,method="spearman"))
with(soma_test,cor.test(mri_severity,soma_mri_severity,method="spearman"))

# Combine training and validation cohorts with predictions
soma_out <- bind_rows(soma_train,soma_test)
