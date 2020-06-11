source("./scripts/pipeline_results_function.R")
source("./scripts/functions.R")
source("./scripts/package-funs.R")

library(gridExtra)
# Interpretation?
library(iml)
# library(ICEbox)

# No adjustments or outliers removed
# p1 <- pipeline_results("msdss_new","MS-DSS Baseline Visit",plat="1.3K",iters=0:119)
# p2 <- pipeline_results("msdss_last_new","MS-DSS Followup Visit",plat="1.3K",iters=0:119)
# jpeg("./results/fig/pipeline_results_20190107.jpeg",units="in",width=12,height=8.5,res=300)
# tiff("./results/fig/pipeline_results_20181109.tiff",units="in",width=12,height
# =8.5,res=300)
# grid.arrange(p1,p2,p3,p4,ncol=2)
# dev.off()

# # Using best MSE
# msdss_new_iter <- 103
# msdss_last_new_iter <- 89

# Using eye ball
msdss_new_iter <- 103
msdss_last_new_iter <- 89

# p1 <- pipeline_results("msdss_new","MS-DSS First Visit (1.3K)",plat="1.3K",iters=0:119,
#                        best_iter = 103)
# p2 <- pipeline_results("msdss_last_new","MS-DSS Last Visit (1.3K)",plat="1.3K",iters=0:119,
#                        best_iter = 89)
# jpeg(my_filename("./results/fig/pipeline_results_eyeball.jpeg"),units="in",width=12,height=8.5,res=300)
# grid.arrange(p1,p2,p3,p4,ncol=2)
# dev.off()

# Age/sex adjustments
# p1 <- pipeline_results("msdss_new_adjust","MS-DSS Baseline Visit - Age/sex Adjusted",plat="1.3K",iters=0:119)
# p2 <- pipeline_results("msdss_last_new_adjust","MS-DSS Followup Visit - Age/sex Adjusted",plat="1.3K",iters=0:119)

# jpeg("./results/fig/adjust_pipeline_results_20190107.jpeg",units="in",width=12,height=8.5,res=300)
# tiff("./results/fig/pipeline_results_20181109.tiff",units="in",width=12,height=8.5,res=300)
# grid.arrange(p1,p2,ncol=2)
# dev.off()

# Using best MSE
# msdss_new_adjust_iter <- 102
# msdss_last_new_adjust_iter <- 106

# p1 <- pipeline_results("msdss_new_adjust","MS-DSS First Visit (1.3K) - Age/sex Adjusted",plat="1.3K",iters=0:119,
#                        best_iter = 97)
# p2 <- pipeline_results("msdss_last_new_adjust","MS-DSS Last Visit (1.3K) - Age/sex Adjusted",plat="1.3K",iters=0:119,
#                        best_iter = 89)
# 
# jpeg(my_filename("./results/fig/adjust_pipeline_results_eyeball.jpeg"),units="in",width=12,height=8.5,res=300)
# grid.arrange(p1,p2,p3,p4,ncol=2)
# dev.off()

# Using eyeball
msdss_new_adjust_iter <- 97
msdss_last_new_adjust_iter <- 89

# p <- pipeline_results("msdss_new_adjust_pms","MS-DSS (P-MS Only)",plat="1.3K",iters=0:119,
#                       best_iter = 102)
msdss_new_adjust_pms_iter <- 102

# p <- pipeline_results("msdss_last_new_adjust_pms","MS-DSS (P-MS Only)",plat="1.3K",iters=0:119,
#                       best_iter = 107)
msdss_last_new_adjust_pms_iter <- 107


# p <- pipeline_results("brain_parfr_new_adjust","Brain Par Fr1",plat="1.3K",iters=0:119,
#                       best_iter = 103)
brain_parfr_new_adjust_iter <- 103

# library(gridExtra)
# Additional Adjustments
# p1 <- pipeline_results("msdss_new_adjust_v2","MS-DSS Baseline Visit",plat="1.3K",iters=0:119)
# p2 <- pipeline_results("msdss_last_new_adjust_v2","MS-DSS Followup Visit",plat="1.3K",iters=0:119)
# p3 <- pipeline_results("brain_parfr_new_adjust_v2","MS-DSS Followup Visit",plat="1.3K",iters=0:119)
# jpeg("./results/fig/pipeline_results_20190107.jpeg",units="in",width=12,height=8.5,res=300)
# tiff("./results/fig/pipeline_results_20181109.tiff",units="in",width=12,height
# =8.5,res=300)
# grid.arrange(p1,p2,p3,p4,ncol=2)
# dev.off()

# # Using best MSE
# msdss_new_iter <- 103
# msdss_last_new_iter <- 89

# Using eye ball
msdss_new_adjust_v2_iter <- 92
msdss_last_new_adjust_v2_iter <- 103
brain_parfr_new_adjust_v2_iter <- 103

# p2 <- pipeline_results("msdss_new_adjust_v2","MS-DSS Baseline Visit",plat="1.3K",iters=0:119,
#                        best_iter=103)
# p2 <- pipeline_results("msdss_last_new_adjust_v2","MS-DSS Followup Visit",plat="1.3K",iters=0:119,
#                        best_iter = 89)
# p3 <- pipeline_results("brain_parfr_new_adjust_v2","MS-DSS Followup Visit",plat="1.3K",iters=0:119,
#                        best_iter = 101)
# jpeg(my_filename("./results/fig/pipeline_results_eyeball.jpeg"),units="in",width=12,height=8.5,res=300)
# grid.arrange(p1,p2,p3,p4,ncol=2)
# dev.off()

iter_list <- c(msdss_new_iter,
               msdss_last_new_iter,
               msdss_new_adjust_iter,
               msdss_last_new_adjust_iter,
               msdss_new_adjust_pms_iter,
               msdss_last_new_adjust_pms_iter,
               brain_parfr_new_adjust_iter,
               msdss_new_adjust_v2_iter,
               msdss_last_new_adjust_v2_iter,
               brain_parfr_new_adjust_v2_iter)
# 
# resp_list <- c("msdss","msdss_last","msdss","msdss_last","msdss","msdss_last","mri_severity","msdss","msdss_last","mri_severity")
# 
# names_list <- c("msdss_new","msdss_last_new","msdss_new_adjust","msdss_last_new_adjust",
#                 "msdss_new_adjust_pms","msdss_last_new_adjust_pms","brain_parfr_new_adjust",
#                 "msdss_new_adjust_v2","msdss_last_new_adjust_v2","brain_parfr_new_adjust")


# iter_list <- c(msdss_new_iter,
#                msdss_last_new_iter,
#                msdss_new_adjust_iter,
#                msdss_last_new_adjust_iter,
#                msdss_new_adjust_pms_iter,
#                msdss_last_new_adjust_pms_iter,
#                brain_parfr_new_adjust_iter)

resp_list <- c("msdss","msdss_last","msdss","msdss_last","msdss","msdss_last","mri_severity",
               "msdss","msdss_last","mri_severity")

names_list <- c("msdss_new","msdss_last_new","msdss_new_adjust","msdss_last_new_adjust",
                "msdss_new_adjust_pms","msdss_last_new_adjust_pms","brain_parfr_new_adjust",
                "msdss_new_adjust_v2","msdss_last_new_adjust_v2","brain_parfr_new_adjust_v2")

fold_list <- c("msdss","msdss_last")

# remove_list <- rep("FALSE",4)

adjust_list <- c(rep("FALSE",2),rep("TRUE",8))

nice_names <- c("MS-DSS First Visit (1.3K)",
                "MS-DSS Last Visit (1.3K)",
                "MS-DSS First Visit (1.3K) - Age/sex Adjusted",
                "MS-DSS Last Visit (1.3K) - Age/sex Adjusted",
                "MS-DSS First Visit (1.3K) - Age/sex Adjusted - PMS",
                "MS-DSS Last Visit (1.3K) - Age/sex Adjusted - PMS",
                "Brain ParFR 1")

platform_list = rep("1.3K",7)

soma_untreated <- read_csv("./data/processed/soma_untreated_20190508.csv")
##################################################################
##################################################################
msdss_dat <- read_csv("Updated_MSDSS_Predictions_10012017_Chris.csv")
train <- filter(soma_untreated,model_cohort=="training")
msdss_train <- filter(msdss_dat,cohort=="training")
sum(msdss_train$patient %in% soma_untreated$patientcode)
sum(msdss_train$patient %in% train$patientcode)
sum(train$patientcode %in% msdss_train$patient)

soma_untreated %>% 
  filter(diagnosis!="HD") %>% 
  group_by(patientcode)  %>% 
  mutate(test = as.numeric(clinic_date_last-min(c(clinic_date,lpdate)))/365.25) %>% 
  .$test %>% 
  mean()
soma_untreated %>% 
  filter(diagnosis!="HD") %>% 
  group_by(patientcode)  %>% 
  mutate(test = as.numeric(clinic_date_last-min(c(clinic_date,lpdate)))/365.25) %>% 
  .$test %>% 
  median()
soma_untreated %>% 
  filter(diagnosis!="HD") %>% 
  group_by(patientcode)  %>% 
  mutate(test = as.numeric(clinic_date_last-min(c(clinic_date,lpdate)))/365.25) %>% 
  .$test %>% 
  sum()



##################################################################
##################################################################

soma_untreated <- soma_untreated %>% 
  mutate(sex = factor(sex,levels=c("Male","Female"))) %>% 
  rename(gender = sex) %>% 
  mutate(gender = factor(ifelse(gender=="Male","M","F"),levels = c("M","F"))) %>% 
  rename(age = age_lp) %>% 
  mutate(age_mri = as.numeric(mri_date-dob)/365)

mri_mod <- lm(volume_calculated_brainparenchymalfr~age_mri,data=soma_untreated)
soma_untreated$mri_severity <- soma_untreated$volume_calculated_brainparenchymalfr - predict(mri_mod,newdata=soma_untreated)

# tiff(my_filename("./results/data_collection/MRI_before.tiff"),res=300,units="in",width=5,height = 4.5)
# tiff(my_filename("./results/data_collection/MRI_before.tiff"),res=300,units="in",width=5.3125,height = 4.25)
# soma_mri <- soma_untreated %>% 
#   filter(!is.na(volume_calculated_brainparenchymalfr))
# 
# mod <- lm(I(1-volume_calculated_brainparenchymalfr)~age_mri,data=soma_mri)
# 
# plot(I(1-volume_calculated_brainparenchymalfr)~age_mri,data=soma_mri,
#      xlab="Age",ylab="Brain Atrophy (1- BPFr)",main="MRI Severity",col=ifelse(residuals(mod)>0,"red","blue"))
# abline(mod)
# dev.off()
#
# # tiff(my_filename("./results/data_collection/MRI_hist.tiff"),res=300,units="in",width=5,height = 4)
# tiff(my_filename("./results/data_collection/MRI_hist.tiff"),res=300,units="in",width=5.3125,height = 4.25)
# hist(I(-1*soma_untreated$mri_severity),nclass=30,col=c(rep("blue",12),rep("red",100)),freq=FALSE,xlab="Atrophy Residuals",main="")
# dev.off()

new <- read_excel("./data/raw/HD Information from Peter - 20181109.xlsx",sheet="1.3k HD")
new <- new[which(!duplicated(with(new,interaction(PatientCode,LPDate)))),]

mod_list <- list()
var_lists <- list()
oob_preds <- list()
new_markers <- str_subset(names(soma_untreated),"^SL")

output <- matrix(0,ncol=length(fold_list),nrow=length(iter_list))
for(i in 1:length(iter_list)){
# for(i in 1:6){
  var_list <- readLines(paste('./scripts/',names_list[i],"/iter",iter_list[i],"_ratios.txt",sep=""))
  var_lists[[i]] <- var_list
}
# saveRDS(var_lists,"./data/processed/all_ratio_list_20190701.rds")
age_to_adjust <- readLines("./data/processed/age_to_explore_20190626.txt")
sex_to_adjust <- readLines("./data/processed/sex_to_explore_20190626.txt")

age_to_adjust_V2 <- readLines("./data/processed/age_to_adjust_V2_20190626.txt")
sex_to_adjust_V2 <- readLines("./data/processed/sex_to_adjust_V2_20190626.txt")

# EDA

marker_list <- lapply(var_lists,function(x){unique(unlist(strsplit(x,"/")))})

sapply(marker_list,function(x){sum(age_to_adjust %in% x)})
sapply(marker_list,function(x){sum(sex_to_adjust %in% x)})

# 
# age_new <- age_to_adjust_V2[age_to_adjust_V2 %notin% age_to_adjust]
# saveRDS(var_lists,my_filename("./data/processed/all_ratio_list.rds"))
new_ratios <- c(var_lists[[1]],var_lists[[2]])

# new_adj_ratios <- c(var_lists[[3]],var_lists[[4]],var_lists[[5]],var_lists[[6]])

new_adj_ratios <- c(var_lists[[3]],var_lists[[4]],var_lists[[5]],var_lists[[6]],
                    var_lists[[7]])

new_adj_v2_ratios <- c(var_lists[[8]],var_lists[[9]],var_lists[[10]])

soma_untreated <- soma_untreated %>% ratio_append(new_ratios,exp=FALSE)

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

soma_train <- soma_untreated %>% filter(model_cohort=="training")
soma_test <- soma_untreated %>% filter(model_cohort=="validation")

soma_adjust_untreated <- soma_untreated %>% ratio_adjust_append(new_adj_ratios,hd_data = new)
soma_adjust_train <- soma_adjust_untreated %>% filter(model_cohort=="training")
soma_adjust_test <- soma_adjust_untreated %>% filter(model_cohort=="validation")


soma_adjust_v2_untreated <- soma_untreated %>% ratio_adjust_append(new_adj_v2_ratios,hd_data = new)
soma_adjust_v2_train <- soma_adjust_v2_untreated %>% filter(model_cohort=="training")
soma_adjust_v2_test <- soma_adjust_v2_untreated %>% filter(model_cohort=="validation")

##############################################################
############### Demographics #################################
##############################################################
new_edss <- read_excel("./data/raw/Patient_EDSS_20190429 - Copy.xlsx")
names(new_edss) <- c("PatientCode","clinic_date","edss_score")
new_edss <- new_edss %>% 
  mutate(clinic_date = ymd(clinic_date)) %>% 
  filter(PatientCode %in% new$PatientCode) %>% 
  remove_dups()
new_msdss <- read_csv("Additional_Endpoints_for_Paav_20190613.csv")
new_msdss <- new_msdss %>% 
  mutate(clinic_date = mdy(clinic_date)) %>% 
  filter(patientcode %in% new$PatientCode) %>% 
  remove_dups() %>% 
  rename(PatientCode = patientcode)

new_msdss <- full_join(new_msdss,new_edss,by=c("PatientCode","clinic_date"))
new_msdss <- new_msdss %>% 
  select(PatientCode,clinic_date,msdss,edss_score)

new_red <- new %>% 
  arrange(PatientCode,LPDate) %>% 
  mutate(LPDate = ymd(LPDate))
new_red <- new_red[!duplicated(new_red$PatientCode),]

new_red <- right_join(new_msdss,new_red,by=c("PatientCode")) %>% 
  mutate(diff = abs(as.numeric(LPDate - clinic_date)/365.25)) %>% 
  group_by(PatientCode,LPDate) %>% 
  filter(diff == min(diff)) %>% 
  ungroup() %>% 
  rename(patientcode = PatientCode,lpdate = LPDate)

soma_all <- bind_rows(soma_train,soma_test,new_red)
soma_all <- soma_all %>% 
  # filter(diagnosis!="HD") %>%
  mutate(dd = as.numeric(lpdate - disease_onset)/365.25)

soma_all %>% 
  group_by(model_cohort) %>% 
  count()

soma_all %>% 
  # filter(diagnosis!="HD") %>% 
  group_by(diagnosis,model_cohort) %>% 
  count() %>% 
  arrange(model_cohort,diagnosis)


soma_all %>% 
  # filter(diagnosis!="HD") %>% 
  group_by(gender,diagnosis,model_cohort) %>% 
  count() %>% 
  arrange(model_cohort,diagnosis,gender)

soma_all %>% 
  # filter(diagnosis!="HD") %>%
  group_by(diagnosis,model_cohort) %>% 
  summarise(mean = mean(age,na.rm=TRUE),sd=sd(age,na.rm=TRUE)) %>% 
  arrange(model_cohort,diagnosis)

soma_all %>% 
  # filter(diagnosis!="HD") %>%
  mutate(dd = as.numeric(lpdate - disease_onset)/365.25) %>%  
  group_by(diagnosis,model_cohort) %>% 
  summarise(mean = mean(dd,na.rm=TRUE),sd=sd(dd,na.rm=TRUE)) %>% 
  arrange(model_cohort,diagnosis)

soma_all %>% 
  # filter(diagnosis!="HD") %>%
  group_by(diagnosis,model_cohort) %>% 
  summarise(mean = mean(edss_score,na.rm=TRUE),sd=sd(edss_score,na.rm=TRUE)) %>% 
  arrange(model_cohort,diagnosis)

soma_all %>% 
  # filter(diagnosis!="HD") %>%
  group_by(diagnosis,model_cohort) %>% 
  summarise(mean = mean(msdss,na.rm=TRUE),sd=sd(msdss,na.rm=TRUE)) %>% 
  arrange(model_cohort,diagnosis)

with(subset(soma_all,diagnosis!="HD"),chisq.test(table(gender,model_cohort)))
with(subset(soma_all,diagnosis!="HD"),chisq.test(table(diagnosis,model_cohort)))
with(subset(soma_all,diagnosis!="HD"),wilcox.test(age~model_cohort))
with(subset(soma_all,diagnosis!="HD"),wilcox.test(dd~model_cohort))
with(subset(soma_all,diagnosis!="HD"),wilcox.test(edss_score~model_cohort))
with(subset(soma_all,diagnosis!="HD"),wilcox.test(msdss~model_cohort))


##############################################################
##############################################################
##############################################################


############################################################
############## Neural Networks #############################
############################################################
# library(keras)
# train_data <- soma_adjust_train %>% 
#   ratio_append(var_lists[[1]]) %>% 
#   select(var_lists[[1]]) %>% 
#   as.matrix()
# train_targets <- soma_adjust_train %>% 
#   .$msdss
# 
# test_data <- soma_adjust_test %>% 
#   ratio_append(var_lists[[1]]) %>% 
#   select(var_lists[[1]]) %>% 
#   as.matrix()
# test_targets <- soma_adjust_test %>% 
#   .$msdss
# 
# mean_train <- apply(train_data, 2, mean)
# std_train <- apply(train_data, 2, sd)
# train_data <- scale(train_data, center = mean_train, scale = std_train)
# test_data <- scale(test_data, center = mean_train, scale = std_train)
# 
# build_model <- function() {
#   model <- keras_model_sequential() %>%
#     layer_dense(units = 32, activation = "relu",
#                 input_shape = dim(train_data)[[2]]) %>%
#     layer_dense(units = 32, activation = "relu") %>%
#     layer_dense(units = 1)
#   
#   model %>% compile(
#     optimizer = "rmsprop",
#     loss = "mse",
#     metrics = c("mae")
#   )
# }
# k <- 4
# indices <- sample(1:nrow(train_data))
# folds <- cut(indices, breaks = k, labels = FALSE)
# 
# num_epochs <- 100
# all_scores <- c()
# for (i in 1:k) {
#   cat("processing fold #", i, "\n")
#   # Prepare the validation data: data from partition # k
#   val_indices <- which(folds == i, arr.ind = TRUE)
#   val_data <- train_data[val_indices,]
#   val_targets <- train_targets[val_indices]
#   
#   # Prepare the training data: data from all other partitions
#   partial_train_data <- train_data[-val_indices,]
#   partial_train_targets <- train_targets[-val_indices]
#   
#   # Build the Keras model (already compiled)
#   model <- build_model()
#   
#   # Train the model (in silent mode, verbose=0)
#   model %>% fit(partial_train_data, partial_train_targets,
#                 epochs = num_epochs, batch_size = 1, verbose = 0)
#   
#   # Evaluate the model on the validation data
#   results <- model %>% evaluate(val_data, val_targets, verbose = 0)
#   all_scores <- c(all_scores, results$mean_absolute_error)
# }
# all_scores
# mean(all_scores)
# # Some memory clean-up
# k_clear_session()
# num_epochs <- 50
# all_mae_histories <- NULL
# for (i in 1:k) {
#   cat("processing fold #", i, "\n")
#   
#   # Prepare the validation data: data from partition # k
#   val_indices <- which(folds == i, arr.ind = TRUE)
#   val_data <- train_data[val_indices,]
#   val_targets <- train_targets[val_indices]
#   
#   # Prepare the training data: data from all other partitions
#   partial_train_data <- train_data[-val_indices,]
#   partial_train_targets <- train_targets[-val_indices]
#   
#   # Build the Keras model (already compiled)
#   model <- build_model()
#   
#   # Train the model (in silent mode, verbose=0)
#   history <- model %>% fit(
#     partial_train_data, partial_train_targets,
#     validation_data = list(val_data, val_targets),
#     epochs = num_epochs, batch_size = 1, verbose = 0
#   )
#   mae_history <- history$metrics$val_mean_absolute_error
#   all_mae_histories <- rbind(all_mae_histories, mae_history)
# }
# average_mae_history <- data.frame(
#   epoch = seq(1:ncol(all_mae_histories)),
#   validation_mae = apply(all_mae_histories, 2, mean)
# )
# library(ggplot2)
# ggplot(average_mae_history, aes(x = epoch, y = validation_mae)) + geom_line()
# ggplot(average_mae_history, aes(x = epoch, y = validation_mae)) + geom_smooth()
# # Get a fresh, compiled model.
# model <- build_model()
# 
# # Train it on the entirety of the data.
# model %>% fit(train_data, train_targets,
#               epochs = 3, batch_size = 16, verbose = 0)
# 
# result <- model %>% evaluate(test_data, test_targets)
# result
# 
# predictions <- model %>% predict(test_data)
# 
# cor.test(test_targets,predictions,method="spearman")
############################################################
############################################################

for(i in 1:length(iter_list)){
# for(i in 1:6){
  if(adjust_list[i]=="FALSE"){soma_mod <- soma_train}
  if(adjust_list[i]=="TRUE" & i %notin% c(8,9,10)){soma_mod <- soma_adjust_train}
  if(adjust_list[i]=="TRUE" & i %in% c(8,9,10)){soma_mod <- soma_adjust_v2_train}
  
  if(i %notin% c(5,6)){
    soma_mod <- soma_mod %>% 
      select(which(names(soma_mod) %in% c(resp_list[i],var_lists[[i]])))
  }
  if(i %in% c(5,6)){
    soma_mod <- soma_mod %>% 
      filter(diagnosis %in% c("PP-MS","SP-MS")) %>% 
      select(which(names(soma_mod) %in% c(resp_list[i],var_lists[[i]])))
  }
  if(i %in% c(7,10)){
    soma_mod <- soma_mod[complete.cases(soma_mod),]
  }
  nvars <- floor(3*sqrt(dim(soma_mod)[2]))
  if(nvars >= dim(soma_mod)[2]-1){nvars <- floor(sqrt(dim(soma_mod)[2]))}
  
  mod <- ranger(data=as.data.frame(soma_mod), num.trees=40000,mtry=nvars,importance='impurity',
                seed=223,dependent.variable.name=resp_list[i],write.forest = TRUE)
  mod_list[[i]] <- mod
}


######################################################################
##################   Add predictions to Dataset ######################
######################################################################
  
soma_train$soma_msdss <- mod_list[[1]]$predictions
soma_train$soma_msdss_last <- mod_list[[2]]$predictions
soma_train$soma_msdss_adjust <- mod_list[[3]]$predictions
soma_train$soma_msdss_last_adjust <- mod_list[[4]]$predictions
# soma_train$soma_msdss_adjust_pms <- mod_list[[5]]$predictions
# soma_train$soma_msdss_last_adjust_pms <- mod_list[[6]]$predictions
# soma_train$soma_mri_severity <- mod_list[[7]]$predictions
soma_train$soma_msdss_adjust_v2 <- mod_list[[8]]$predictions
soma_train$soma_msdss_last_adjust_v2 <- mod_list[[9]]$predictions

soma_train_mri <- filter(soma_train,!is.na(mri_severity))
soma_train_mri$soma_mri_severity <- (-1)*mod_list[[7]]$predictions
soma_train_mri <- select(soma_train_mri,sampleid,soma_mri_severity)
soma_train <- left_join(soma_train,soma_train_mri,by=c("sampleid"))


soma_train_mri <- filter(soma_train,!is.na(mri_severity))
soma_train_mri$soma_mri_severity_v2 <- (-1)*mod_list[[10]]$predictions
soma_train_mri <- select(soma_train_mri,sampleid,soma_mri_severity_v2)
soma_train <- left_join(soma_train,soma_train_mri,by=c("sampleid"))


soma_test$soma_msdss <- predict(mod_list[[1]],soma_test)$predictions
soma_test$soma_msdss_last <- predict(mod_list[[2]],soma_test)$predictions
soma_test$soma_msdss_adjust <- predict(mod_list[[3]],soma_adjust_test)$predictions
soma_test$soma_msdss_last_adjust <- predict(mod_list[[4]],soma_adjust_test)$predictions
soma_test$soma_msdss_adjust_pms <- predict(mod_list[[5]],soma_adjust_test)$predictions
soma_test$soma_msdss_last_adjust_pms <- predict(mod_list[[6]],soma_adjust_test)$predictions
soma_test$soma_mri_severity <- predict(mod_list[[7]],soma_adjust_test)$predictions

soma_test$soma_msdss_adjust_v2 <- predict(mod_list[[8]],soma_adjust_v2_test)$predictions
soma_test$soma_msdss_last_adjust_v2 <- predict(mod_list[[9]],soma_adjust_v2_test)$predictions
soma_test$soma_mri_severity_v2 <- predict(mod_list[[10]],soma_adjust_v2_test)$predictions


soma_test$soma_mri_severity <- (-1)*soma_test$soma_mri_severity
soma_test$soma_mri_severity_v2 <- (-1)*soma_test$soma_mri_severity_v2
soma_test$mri_severity <- (-1)*soma_test$mri_severity
soma_train$mri_severity <- (-1)*soma_train$mri_severity

pred_list <- c("soma_msdss_adjust_v2","soma_msdss_last_adjust_v2","soma_mri_severity_v2")

# out_dat <- bind_rows(soma_train,soma_test) %>%
#   select(sampleid,patientcode,pred_list)
# write_csv(out_dat,"Model Predictions for Heatmap - 20190701.csv")

library(DescTools)
my_CCC<- function(x,y,digits=3){
  # digits <- 3
  obj <- CCC(x,y)
  obj <- obj[[1]]
  est <- round(obj[1],digits)
  li <- round(obj[2],digits)
  ui <- round(obj[3],digits)
  out <- paste(est," (",li," - ",ui,")",sep="")
  return(out)
}

############################################################################
############### Variable Importance/ Ratio correlations for 1.3K Models ####
############################################################################
# invisible(lapply(1:7,function(i){
#   mod <- mod_list[[i]]
#   imp_dat <- cbind(clean_rat(names(mod$variable.importance)),importance=as.numeric(mod$variable.importance))
#   imp_dat <- imp_dat %>%
#     mutate(msdss_rho = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_untreated[[xchar]],soma_adjust_untreated$msdss,method="spearman")$est)})) %>%
#     mutate(msdss_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_untreated[[xchar]],soma_adjust_untreated$msdss,method="spearman")$p.value)})) %>%
#     mutate(msdss_pvalue = p.adjust(msdss_pvalue,method="fdr")) %>%
#     mutate(mri_severity_rho = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_untreated[[xchar]],soma_adjust_untreated$mri_severity,method="spearman")$est)})) %>%
#     mutate(mri_severity_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_untreated[[xchar]],soma_adjust_untreated$mri_severity,method="spearman")$p.value)})) %>%
#     mutate(mri_severity_pvalue = p.adjust(mri_severity_pvalue,method="fdr")) %>%
#     write_csv(my_filename(paste("./results/importance_added/",names_list[i],"_ratios.csv",sep="")))
# }))

# invisible(lapply(8:10,function(i){
#   mod <- mod_list[[i]]
#   imp_dat <- cbind(clean_rat(names(mod$variable.importance)),importance=as.numeric(mod$variable.importance))
#   imp_dat <- imp_dat %>%
#     mutate(msdss_rho = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$est)})) %>%
#     mutate(msdss_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$p.value)})) %>%
#     mutate(msdss_pvalue = p.adjust(msdss_pvalue,method="fdr")) %>%
#     mutate(msdss_last_rho = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$est)})) %>%
#     mutate(msdss_last_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$p.value)})) %>%
#     mutate(msdss_last_pvalue = p.adjust(msdss_last_pvalue,method="fdr")) %>%
#     mutate(mri_severity_rho = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$est)})) %>%
#     mutate(mri_severity_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$p.value)})) %>%
#     mutate(mri_severity_pvalue = p.adjust(mri_severity_pvalue,method="fdr")) %>%
#     write_csv(my_filename(paste("./results/importance_added/",names_list[i],"_ratios.csv",sep="")))
# }))


# marker_clusters <- read_excel("./results/string/pipeline_clusters/Pipeline_Markers_Clusters_20190702.xlsx")
# invisible(lapply(8:10,function(i){
#   # i <- 8
#   mod <- mod_list[[i]]
#   imp_dat <- cbind(clean_rat(names(mod$variable.importance)),importance=as.numeric(mod$variable.importance))
#   imp_dat <- imp_dat %>%
#     mutate(importance = 100*importance/sum(importance)) %>% 
#     left_join(select(marker_clusters,v1=marker,cluster1 = cluster),by="v1") %>% 
#     select(ratio,v1,gene1,protein1,cluster1,everything()) %>% 
#     left_join(select(marker_clusters,v2=marker,cluster2 = cluster),by="v2") %>% 
#     select(ratio,v1,gene1,protein1,cluster1,v2,gene2,protein2,cluster2,everything()) %>% 
#     
#     mutate(msdss_rho = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$est)})) %>%
#     mutate(msdss_v1_rho = sapply(v1,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$est)})) %>%
#     mutate(msdss_v2_rho = sapply(v2,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$est)})) %>%
#     mutate(msdss_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$p.value)})) %>%
#     mutate(msdss_pvalue = p.adjust(msdss_pvalue,method="fdr")) %>%
#     mutate(msdss_v1_pvalue = sapply(v1,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$p.value)})) %>%
#     mutate(msdss_v1_pvalue = p.adjust(msdss_v1_pvalue,method="fdr")) %>%
#     mutate(msdss_v2_pvalue = sapply(v2,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss,method="spearman")$p.value)})) %>%
#     mutate(msdss_v2_pvalue = p.adjust(msdss_v2_pvalue,method="fdr")) %>%
#     
#     mutate(msdss_last_rho = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$est)})) %>%
#     mutate(msdss_last_v1_rho = sapply(v1,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$est)})) %>%
#     mutate(msdss_last_v2_rho = sapply(v2,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$est)})) %>%
#     mutate(msdss_last_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$p.value)})) %>%
#     mutate(msdss_last_pvalue = p.adjust(msdss_last_pvalue,method="fdr")) %>%
#     mutate(msdss_last_v1_pvalue = sapply(v1,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$p.value)})) %>%
#     mutate(msdss_last_v1_pvalue = p.adjust(msdss_last_v1_pvalue,method="fdr")) %>%
#     mutate(msdss_last_v2_pvalue = sapply(v2,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$msdss_last,method="spearman")$p.value)})) %>%
#     mutate(msdss_last_v2_pvalue = p.adjust(msdss_last_v2_pvalue,method="fdr")) %>%
#     
#     mutate(mri_severity_rho = (-1)*sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$est)})) %>%
#     mutate(mri_severity_v1_rho = (-1)*sapply(v1,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$est)})) %>%
#     mutate(mri_severity_v2_rho = (-1)*sapply(v2,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$est)})) %>%
#     mutate(mri_severity_pvalue = sapply(ratio,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$p.value)})) %>%
#     mutate(mri_severity_pvalue = p.adjust(mri_severity_pvalue,method="fdr")) %>%
#     mutate(mri_severity_v1_pvalue = sapply(v1,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$p.value)})) %>%
#     mutate(mri_severity_v1_pvalue = p.adjust(mri_severity_v1_pvalue,method="fdr")) %>%
#     mutate(mri_severity_v2_pvalue = sapply(v2,function(xchar){as.numeric(cor.test(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated$mri_severity,method="spearman")$p.value)})) %>%
#     mutate(mri_severity_v2_pvalue = p.adjust(mri_severity_v2_pvalue,method="fdr")) %>%
#     
#     write_csv(my_filename(paste("./results/importance_added/",names_list[i],"_ratios.csv",sep="")))
# }))

################################################################################
######################### Grouped by procedure #################################
################################################################################
my_scatterplot <- function(x,y,main="",fix=TRUE,xlab=NULL,ylab=NULL,cor_method="spearman",reg=FALSE,
                           xloc = 4, yloc = 1, ...){
  if(fix==TRUE){
    if(is.null(xlab)){xlab <- "Predicted"}
    if(is.null(ylab)){ylab <- "Observed"}
    lim <- range(x,y)
    plot(x,y,xlim=lim,ylim=lim,xlab=xlab,ylab=ylab,
         main=main,cex.lab=2,cex.main=2,cex.axis=1.3,...)
    # text(4,1.5,my_CCC(x,y),cex=2)
    # text(xloc,yloc,paste("Spearman =",my_cor(x,y,method=cor_method)),cex=2)
    # legend("bottomright",paste("Rho =",my_cor(x,y,method=cor_method)),cex=2,bty="n")
    if(reg==TRUE){
      abline(lm(y~x))
    }else{abline(0,1)}
  }  
  if(fix==FALSE){
    plot(x,y,xlab=xlab,ylab=ylab,
         main=main,cex.lab=2,cex.main=2,...)
    # text(xloc,yloc,paste("Rho =",my_cor(x,y,method=cor_method)),cex=2)
    # legend("topleft",paste("Spearman =",my_cor(x,y,method=cor_method)),cex=2,bty="n")
    if(reg==TRUE){
      abline(lm(y~x))
    }else{abline(0,1)}}
}

# jpeg(my_filename("./results/fig/Validation_Predictions.JPG"),units="in",res=300,width=8.5,height=12)
# par(mfrow=c(4,2))
# with(soma_test,my_scatterplot(soma_msdss,msdss,main="MS-DSS Baseline Visit"))
# with(soma_test,my_scatterplot(soma_msdss_last,msdss_last,main="MS-DSS Followup Visit"))
# with(soma_test,my_scatterplot(soma_msdss_adjust,msdss,main="MS-DSS Baseline Visit - Adjusted"))
# with(soma_test,my_scatterplot(soma_msdss_last_adjust,msdss_last,main="MS-DSS Followup Visit - Adjusted"))
# with(subset(soma_test,diagnosis %in% c("PP-MS","SP-MS")),
#      my_scatterplot(soma_msdss_adjust_pms,msdss,main="MS-DSS Baseline Visit - Adjusted - PMS"))
# with(subset(soma_test,diagnosis %in% c("PP-MS","SP-MS")),
#      my_scatterplot(soma_msdss_last_adjust_pms,msdss_last,main="MS-DSS Followup Visit - Adjusted - PMS"))
# with(subset(soma_test,!is.na(mri_severity)),my_scatterplot(soma_mri_severity,mri_severity,main="MRI Severity - Adjusted"))
# par(mfrow=c(1,1))
# dev.off()
# 
# par(mfrow=c(2,3))
# with(soma_test,my_scatterplot(soma_msdss_adjust,msdss,main="MS-DSS Baseline Visit - Initial"))
# with(soma_test,my_scatterplot(soma_msdss_last_adjust,msdss_last,main="MS-DSS Followup Visit - Initial"))
# with(subset(soma_test,!is.na(mri_severity)),my_scatterplot(soma_mri_severity,mri_severity,main="MRI Severity - Initial"))
# 
# with(soma_test,my_scatterplot(soma_msdss_adjust_v2,msdss,main="MS-DSS Baseline Visit - Updated"))
# with(soma_test,my_scatterplot(soma_msdss_last_adjust_v2,msdss_last,main="MS-DSS Followup Visit - Updated"))
# with(subset(soma_test,!is.na(mri_severity)),my_scatterplot(soma_mri_severity_v2,mri_severity,main="MRI Severity - Updated"))
# par(mfrow=c(1,1))

# jpeg(my_filename("./results/fig/Validation_Predictions_Final.JPG"),units="in",res=300,width=8.5,height=6)
tiff(my_filename("./results/pipeline_results/Validation_Predictions_Final.tiff"),units="in",res=300,width=13,height=5)
par(mfrow=c(1,3),mar = c(2, 2, 2, 2) + 0.1)
with(soma_test,my_scatterplot(soma_msdss_adjust_v2,msdss,main="MS-DSS Baseline",pch=16))
with(soma_test,my_scatterplot(soma_msdss_last_adjust_v2,msdss_last,main="MS-DSS Follow-up",pch=16))
with(subset(soma_test,!is.na(mri_severity)),my_scatterplot(soma_mri_severity_v2,mri_severity,main="MRI Severity",pch=16))
par(mfrow=c(1,1),mar = c(5, 4, 4, 2) + 0.1)
dev.off()

# jpeg(my_filename("./results/fig/Validation_Predictions_Final.JPG"),units="in",res=300,width=8.5,height=6)
# tiff(my_filename("./results/pipeline_results/Validation_Predictions_Final.tiff"),units="in",res=300,width=13,height=5)
tiff(my_filename("./results/pipeline_results/Validation_Predictions_Final.tiff"),units="in",res=300,width=7.8,height=3)
par(mfrow=c(1,3),mar = c(2, 2, 2, 2) + 0.1)
with(soma_test,my_scatterplot(soma_msdss_adjust_v2,msdss,main=""))
with(subset(soma_test,!is.na(mri_severity)),my_scatterplot(soma_mri_severity_v2,mri_severity,main=""))
with(soma_test,my_scatterplot(soma_msdss_last_adjust_v2,msdss_last,main=""))
par(mfrow=c(1,1),mar = c(5, 4, 4, 2) + 0.1)
dev.off()

tiff(my_filename("./results/pipeline_results/Validation_Predictions_Baseline.tiff"),units="in",res=300,width=4.5,height=5)
with(soma_test,my_scatterplot(soma_msdss_adjust_v2,msdss,main="",
                              xlab="",ylab=""))
dev.off()

tiff(my_filename("./results/pipeline_results/Validation_Predictions_Followup.tiff"),units="in",res=300,width=4.5,height=5)
with(soma_test,my_scatterplot(soma_msdss_last_adjust_v2,msdss_last,main="",
                              xlab="",ylab=""))
dev.off()

tiff(my_filename("./results/pipeline_results/Validation_Predictions_MRI.tiff"),units="in",res=300,width=4.5,height=5)
with(subset(soma_test,!is.na(mri_severity)),my_scatterplot(soma_mri_severity_v2,mri_severity,main="",
                                                           xlab="",ylab=""))
dev.off()

############################################################################
######################### Outcome Data #####################################
############################################################################
marker_clusters <- read_excel("./results/string/pipeline_clusters/Pipeline_Markers_Clusters_20190702.xlsx")

new_outcome <- bind_rows(soma_test,soma_train)
new_outcome <- new_outcome %>% 
  select(sampleid,starts_with("soma")) %>% 
  select(-soma_id)

soma_untreated <- left_join(soma_untreated,new_outcome,by="sampleid")

outcomes <- c("msss","armss","msdss")
outcomes <- c(outcomes,paste(outcomes,"last",sep="_"))
# outcomes <- c(outcomes,c("combiwise_score_slope_unadj","combiwise_score_slope_adj","neurex_total_slope_unadj","neurex_total_slope_adj"))
outcomes <- c(outcomes,c("combiwise_score_slope_adj","neurex_total_slope_adj"))

# outcomes <- outcomes[-c(10)]
# volumes <- str_subset(names(soma_untreated),"volume")
# volumes <- names(which(apply(soma_untreated[,volumes],2,function(x){sum(is.na(x))})==23))
# volumes <- volumes[!str_detect(volumes,"source")]
# outcomes <- c(outcomes,volumes,"mri_severity")

soma_markers <- names(new_outcome)[-1]
# soma_markers <- soma_markers[1:4]

soma_selected <- soma_untreated %>% 
  filter(model_cohort=="validation")

# 1.3K Platform, Age/sex adjusted somamers
cor_est <- function(x,y,...){
  obj <- cor.test(x,y,...)
  val <- obj$estimate
  return(as.numeric(val))
}
cor_pval <- function(x,y,...){
  obj <- cor.test(x,y,...)
  pval <- obj$p.value
  return(as.numeric(pval))
}

##########################################################
################### For Erin #############################
##########################################################
module1_markers <- marker_clusters$marker[marker_clusters$cluster==1]
all_rats <- unique(unlist(var_lists[8:10]))


module1_rats <- all_rats[sapply(all_rats,function(x){any(str_detect(x,module1_markers))})]
module1 <- c(module1_markers,module1_rats)

module1_est <- sapply(module1,function(ychar){
  sapply(c(outcomes,"mri_severity"),function(xchar){
    cor_est(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated[[ychar]],method="spearman")
  })
})
module1_pval <- sapply(module1,function(ychar){sapply(c(outcomes,"mri_severity"),function(xchar){
  cor_pval(soma_adjust_v2_untreated[[xchar]],soma_adjust_v2_untreated[[ychar]],method="spearman")})})

module1_est <- t(module1_est)
module1_pval <- t(module1_pval)

module1_pval <- as_tibble(module1_pval) %>% 
  mutate_all(function(x){p.adjust(x,method="fdr")}) %>%
  # mutate(outcome = c(outcomes,"mri_severity"),type="p-value") %>% 
  # select(outcome,type,everything())
  mutate(marker=module1,type="p-value") %>% 
  select(marker,type,everything())

module1_est <- as_tibble(module1_est) %>% 
  # mutate(outcome = c(outcomes,"mri_severity"),type="spearman") %>% 
  # select(outcome,type,everything())
  mutate(marker = module1,type="spearman") %>% 
  select(marker,type,everything())

module1_dat <- bind_rows(module1_est,module1_pval)
module1_dat <- module1_dat %>% 
  # gather(key,value,module1) 
  gather(key,value,outcomes,"mri_severity") 

module1_dat <- module1_dat %>% 
  mutate(key = factor(key,levels=unique(module1_dat$key))) %>% 
  mutate(type=factor(type,levels=c("spearman","p-value"))) %>% 
  # arrange(outcome,key,type) %>% 
  arrange(marker,key,type) %>% 
  mutate(key = paste(as.character(key),type)) %>% 
  select(-type) 

module1_dat <- module1_dat %>% 
  mutate(key = factor(key,levels=unique(module1_dat$key))) %>% 
  spread(key,value)

module1_marker_dat <- module1_dat %>% 
  filter(nchar(marker)==8)
# module1_marker_dat <- as_tibble(cbind(clean_marker(module1_marker_dat[[1]]),module1_marker_dat[,-1]))
module1_marker_dat <- right_join(marker_clusters,module1_marker_dat,by="marker")

module1_ratio_dat <- module1_dat %>% 
  filter(nchar(marker)==17)
module1_ratio_dat <- as_tibble(cbind(clean_rat(module1_ratio_dat[[1]]),module1_ratio_dat[,-1]))
names(module1_ratio_dat)[c(1:5)] <- c("ratio","marker","gene","protein","marker2")

module1_ratio_dat <- left_join(module1_ratio_dat,select(marker_clusters,marker,cluster,nice_name),by="marker")
module1_ratio_dat <- left_join(module1_ratio_dat,select(marker_clusters,marker2=marker,cluster2=cluster,nice_name2=nice_name),by="marker2")
module1_ratio_dat <- module1_ratio_dat %>% 
  select(ratio,marker,gene,protein,cluster,nice_name,marker2,gene2,protein2,cluster2,nice_name2,everything())

module1_dat <- bind_rows(module1_ratio_dat,module1_marker_dat)
# write_csv(module1_dat,"Module1 Data for Erin_20190725.csv",na = "")


##########################################################
##########################################################

new_adj_est <- sapply(soma_markers,function(ychar){
  sapply(outcomes,function(xchar){
    cor_est(soma_selected[[xchar]],soma_selected[[ychar]],method="spearman")
  })
})
new_adj_pval <- sapply(soma_markers,function(ychar){sapply(outcomes,function(xchar){
  cor_pval(soma_selected[[xchar]],soma_selected[[ychar]],method="spearman")})})



new_adj_pval <- as_tibble(new_adj_pval) %>% 
  mutate_all(function(x){p.adjust(x,method="fdr")}) %>%
  mutate(outcome = outcomes,type="p-value") %>% 
  select(outcome,type,everything())

new_adj_est <- as_tibble(new_adj_est) %>% 
  mutate(outcome = outcomes,type="spearman") %>% 
  select(outcome,type,everything())

new_adj <- bind_rows(new_adj_est,new_adj_pval)
new_adj <- new_adj %>% 
  gather(key,value,soma_markers) 

new_adj <- new_adj %>% 
  mutate(key = factor(key,levels=unique(new_adj$key))) %>% 
  mutate(type=factor(type,levels=c("spearman","p-value"))) %>% 
  arrange(outcome,key,type) %>% 
  mutate(key = paste(as.character(key),type)) %>% 
  select(-type) 

new_adj <- new_adj %>% 
  mutate(key = factor(key,levels=unique(new_adj$key))) %>% 
  spread(key,value)

# write_csv(new_adj,my_filename("Somalogic Severity Model vs Severity Outcomes.csv"))
# write_csv(new_adj,my_filename("./results/importance_added/Somalogic Severity Model vs Severity Outcomes.csv"))

# par(mfrow=c(length(outcomes),3))
# for(i in outcomes){
#   for(j in soma_markers[8:10]){
#     plot(soma_selected[[j]],soma_selected[[i]],xlab=j,ylab=i,pch=16) 
#   }
# }
# par(mfrow=c(1,1))
# outcomes_to_plot <- outcomes[-c(1,5,8,9,11)]
# par(mfrow=c(3,length(outcomes_to_plot)),mar = c(5, 4, 4, 2) + 0.1)
# # par(mar=c(0,0,0,0) + 0.1)
# for(i in soma_markers[8:10]){
#   for(j in outcomes_to_plot){
#     plot(soma_selected[[j]],soma_selected[[i]],xlab="",ylab="",pch=16,xaxt="n",yaxt="n")
#   }
# }
# par(mfrow=c(1,1),mar = c(5, 4, 4, 2) + 0.1)

corr_out <- matrix(0,nrow=length(outcomes_to_plot),ncol = 3)
colnames(corr_out) <- soma_markers[8:10]
rownames(corr_out) <- outcomes_to_plot
for(i in soma_markers[8:10]){
  for(j in outcomes_to_plot){
    corr_out[j,i] <- as.numeric(cor.test(soma_selected[[j]],soma_selected[[i]],method="spearman")$est)
  }
}

library(ComplexHeatmap)
library(circlize)
Heatmap(corr_out,col=colorRamp2(c(-1,0,1),c("blue","white","red"),space="sRGB"),cluster_rows = FALSE,cluster_columns = FALSE)

cor_col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061")))(200)
corrplot(corr_out,col = cor_col,sig.level = .05,stars=TRUE)


add_lines <- function(x,y,line_col = "red"){
  mod <- lm(y~x)
  # abline(mod)
  xind <- seq(from = min(x),to = max(x),length=1000)
  upper <- predict(mod,newdata = data.frame(x=xind),interval="confidence")[,3]
  pred <- predict(mod,newdata = data.frame(x=xind),interval="confidence")[,1]
  lower <- predict(mod,newdata = data.frame(x=xind),interval="confidence")[,2]
  for(i in 1:length(xind)){polygon(rep(xind[i],2),c(upper[i],lower[i]),density=NA,col="grey")}
  lines(xind,pred,lty=1,lwd=2,col=line_col)
  # lines(xind,upper,lty=2,lwd=2,col=2)
  # lines(xind,lower,lty=2,lwd=2,col=2)
}

tiff(my_filename("./results/pipeline_results/correlations_with_outcomes.tiff"),res=300,units="in",width=8.5,height=11)
# outcomes_to_plot <- c("msss","msss_last","armss","armss_last","msdss","combiwise_score_slope_adj","neurex_total_slope_adj")
outcomes_to_plot <- c("msss","msss_last","armss","armss_last","msdss","msdss_last","neurex_total_slope_adj","combiwise_score_slope_adj")
par(mfrow=c(length(outcomes_to_plot),3),mar = c(5, 4, 4, 2) + 0.1)
par(mar=c(0,0,0,0) + 2)
for(i in  outcomes_to_plot){
  for(j in soma_markers[c(8,10,9)]){
    est <- round(as.numeric(new_adj[new_adj$outcome==i,paste(j,"spearman")]),3)
    pval <- round(as.numeric(new_adj[new_adj$outcome==i,paste(j,"p-value")]),3)
    main_col <- ifelse(pval<=0.05,"red","black")
    pval <- format.pval(pval,eps=0.001)
    if(j != "soma_msdss_adjust_v2" & i != "combiwise_score_slope_adj"){
      plot(soma_selected[[j]],soma_selected[[i]],xlab="",ylab="",xaxt="n",yaxt="n",main=paste("Rho = ",est," (",pval,")",sep=""),
           col.main=main_col,tck=0.05,type="n")
      axis(side=1,labels=FALSE,tck=0.05)
      axis(side=2,labels=FALSE,tck=0.05)
      add_lines(soma_selected[[j]],soma_selected[[i]],line_col = main_col)
      points(soma_selected[[j]],soma_selected[[i]])
    }
    if(j == "soma_msdss_adjust_v2" & i == "combiwise_score_slope_adj"){
      plot(soma_selected[[j]],soma_selected[[i]],xlab="",ylab="",main=paste("Rho = ",est," (",pval,")",sep=""),
           col.main=main_col,tck=0.05,type="n")
      axis(side=1,labels=FALSE,tck=0.05)
      axis(side=2,labels=FALSE,tck=0.05)
      add_lines(soma_selected[[j]],soma_selected[[i]],line_col = main_col)
      points(soma_selected[[j]],soma_selected[[i]])
    }
    if(j == "soma_msdss_adjust_v2" & i != "combiwise_score_slope_adj"){
      plot(soma_selected[[j]],soma_selected[[i]],xlab="",ylab="",xaxt="n",main=paste("Rho = ",est," (",pval,")",sep=""),
           col.main=main_col,tck=0.05,type="n")
      add_lines(soma_selected[[j]],soma_selected[[i]],line_col = main_col)
      axis(side=1,labels=FALSE,tck=0.05)
      axis(side=2,labels=FALSE,tck=0.05)
      points(soma_selected[[j]],soma_selected[[i]])
    }
    if(j != "soma_msdss_adjust_v2" & i == "combiwise_score_slope_adj"){
      plot(soma_selected[[j]],soma_selected[[i]],xlab="",ylab="",yaxt="n",main=paste("Rho = ",est," (",pval,")",sep=""),
           col.main=main_col,tck=0.05,type="n")
      axis(side=1,labels=FALSE,tck=0.05)
      axis(side=2,labels=FALSE,tck=0.05)
      add_lines(soma_selected[[j]],soma_selected[[i]],line_col = main_col)
      points(soma_selected[[j]],soma_selected[[i]])
    }
  }
}
par(mfrow=c(1,1),mar = c(5, 4, 4, 2) + 0.1)
dev.off()
