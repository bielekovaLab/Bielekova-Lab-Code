library(ggplot2)
library(multcomp)
library(effects)

source("./scripts/functions.R")
source("./scripts/package-funs.R")

library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(dendextend)
library(dynamicTreeCut)
library(colorspace)
library(RColorBrewer)
library(gridExtra)
library(ggbeeswarm)

z_score <- function(x,log=FALSE){
  if(log==TRUE){x <- log(x)}
  out <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm = TRUE)
  return(out)
}

var_lists <- readRDS("./data/processed/all_ratio_list_20190701.rds")
# var_lists <- readRDS("./data/processed/all_ratio_list_20190515.rds")
# age_to_adjust <- readLines("./data/processed/age_to_explore_20190510.txt")
# sex_to_adjust <- readLines("./data/processed/sex_to_explore_20190510.txt")
# age_to_adjust <- readLines("./data/processed/age_to_explore_20190626.txt")
# sex_to_adjust <- readLines("./data/processed/sex_to_explore_20190626.txt")
age_to_adjust <- readLines("./data/processed/age_to_adjust_V2_20190626.txt")
sex_to_adjust <- readLines("./data/processed/sex_to_adjust_V2_20190626.txt")

############################
# sum(age_to_adjust %in% unique(unlist(strsplit(unlist(var_lists[c(3:4,7)]),"/"))))

############################

model_predictions <- read_csv("Model Predictions for Heatmap - 20190701.csv")

soma_untreated <- read_csv("./data/processed/soma_untreated_20190508.csv") 
soma_untreated <- soma_untreated %>% 
  rename(age=age_lp) %>% 
  mutate(age_mri = as.numeric(mri_date-dob)/365)

soma_untreated$armss_last[soma_untreated$patientcode=="NIB371"] <- 10
soma_untreated$armss_last[soma_untreated$patientcode=="NIB029"] <- 2.11

soma_untreated <- soma_untreated %>% 
  mutate(volume_calculated_lesionfr = volume_lesions/volume_calculated_totalbrain) %>% 
  mutate(volume_calculated_csf_fr = volume_calculated_totalcsf/volume_calculated_totalbrain)

# soma_untreated <- soma_untreated %>% 
#   mutate(atr = 1-volume_calculated_brainparenchymalfr)

summary(mri_mod <- lm(volume_calculated_brainparenchymalfr~age_mri,data=soma_untreated))
# summary(mri_mod_check <- lm(atr~age_mri,data=soma_untreated))

soma_untreated$mri_severity <- (-1)*(soma_untreated$volume_calculated_brainparenchymalfr - predict(mri_mod,newdata=soma_untreated))
# soma_untreated$mri_severity_check <- soma_untreated$atr - predict(mri_mod_check,newdata=soma_untreated)

diffs  <- c("t1_black_hole_fraction","t2_lesion_load",
            "atrophy_type","cel_lesion_load","patient","mri_date","cel_exact_number",
            "atrophy_score","atrophy_score_divide_by_cerebrum",
            "lesion_score","lesion_score_plus_cerebrum","gm_plus_juxtacortical")
clean_levels <- function(x,factorize = TRUE){
  out <- str_replace_all(unlist(lapply(strsplit(x,"[(]"),first))," ","")
  if(factorize==TRUE){out <- ordered(out,levels=c("None","Mild","Moderate","Severe"))}
  return(out)
}

# qual <- c("t1_black_hole_fraction","t2_lesion_load","cel_lesion_load","juxtacortical_lesion_load",
#           "lesion_atrophy_brainstem","lesion_load_brainstem","gm_lesion_load",
#           "lesion_load_cerebellum","lesion_atrophy_cerebellum",
#           "lesion_atrophy_medulla",
#           "lesion_load_medulla",
#           "lesion_cerebrum")
qual <- c("t1_black_hole_fraction","t2_lesion_load","cel_lesion_load","juxtacortical_lesion_load",
          "lesion_atrophy_brainstem","lesion_load_brainstem","gm_lesion_load",
          "lesion_load_cerebellum","lesion_atrophy_cerebellum",
          "lesion_atrophy_medulla",
          "lesion_load_medulla")

soma_untreated <- soma_untreated %>% 
  mutate_at(names(soma_untreated)[names(soma_untreated) %in% qual& names(soma_untreated) %notin% diffs],clean_levels)
soma_untreated$t1_black_hole_fraction <- ordered(soma_untreated$t1_black_hole_fraction,levels=c("None","<50%",">50%")) 
soma_untreated$t2_lesion_load <- ordered(clean_levels(soma_untreated$t2_lesion_load,factorize = FALSE),
                                        levels = c("None","Mild","Moderate","Severe"))
# soma_untreated$atrophy_type <- factor(soma_untreated$atrophy_type)
soma_untreated$cel_lesion_load <- ordered(clean_levels(soma_untreated$cel_lesion_load,factorize = FALSE),
                                         levels=c("None","Mild","Moderate","Severe"))
soma_untreated <- soma_untreated  %>% 
  mutate(csf_albuminquotient = csf_albuminquotient - 3.0496 + 0.0493*age)


soma_new <- soma_untreated %>% 
  rename(gender=sex) %>% 
  mutate(gender = ifelse(gender=="Male","M","F")) 

new_markers <- soma_untreated %>% 
  select(starts_with("SL")) %>% 
  names()

new <- read_excel("./data/raw/HD Information from Peter - 20181109.xlsx",sheet="1.3k HD")
new <- new[which(!duplicated(with(new,interaction(PatientCode,LPDate)))),]
names(new)[1:10] <- tolower(names(new))[1:10]
new <- new %>% 
  mutate(gender = factor(gender,levels=c("M","F"))) 

soma_untreated_adj <- soma_untreated %>%
  select(-starts_with("SL"))

# Add adjusted proteins to dataset
soma_new_adj <- foreach(i=1:length(new_markers),.final=as.data.table,.packages=c("dplyr")) %do% {
  mod <- adjust_function(new_markers[i],data=new)
  (log(soma_new[[new_markers[i]]]) - predict(mod,newdata=soma_new))
}
names(soma_new_adj) <- new_markers
soma_new_adj[,sampleid := soma_new$sampleid]
setcolorder(soma_new_adj,c("sampleid"))
soma_untreated_adj <- merge(soma_untreated_adj,soma_new_adj,by="sampleid")


new_adj <- foreach(i=1:length(new_markers),.final=as.data.table,.packages=c("dplyr")) %do% {
  mod <- adjust_function(new_markers[i],data=new)
  (log(new[[new_markers[i]]]) - predict(mod,newdata=new))
}
names(new_adj) <- new_markers
new_adj[,sampleid := new$sampleid]
setcolorder(new_adj,c("sampleid"))
# 


ind <- c(3:10)
pipeline_rats <- unique(unlist(var_lists[ind]))

soma_untreated <- soma_untreated %>% 
  ratio_append(pipeline_rats,adj=FALSE)

new <- new %>% 
  mutate(lpdate = ymd(lpdate)) %>% 
  ratio_append(pipeline_rats,adj=FALSE)

extra_csf <- clean(read_excel("./data/raw/Patient_NIHLabCSF1_20190508.xlsx"))
extra_csf <- extra_csf %>% 
  filter(!is.na(patientcode)&!is.na(protocol1)&!is.na(date)) %>% 
  select(-protocol1) %>% 
  mutate(date = ymd(date)) %>% 
  rename(lpdate = date) %>% 
  remove_dups() %>% 
  mutate(id = paste(patientcode,lpdate)) 

new <- right_join(extra_csf,new,by=c("patientcode","lpdate"))
new$csf_ocb <- as.integer(new$csf_ocb) # Dumb.....

soma_untreated <- soma_untreated %>% 
  bind_rows(new) 

soma_untreated <- soma_untreated %>% 
  mutate_at(vars(starts_with("SL")),z_score,log=TRUE)

new <- new %>% 
  select(-starts_with("SL")) %>% 
  inner_join(new_adj,by="sampleid")

soma_untreated_adj <- soma_untreated_adj %>% 
  bind_rows(new) %>% 
  ratio_append(pipeline_rats,adj = TRUE) %>% 
  mutate_at(vars(starts_with("SL")),z_score,log=FALSE)

soma_untreated <- soma_untreated %>% 
  mutate(ms_diagnosis = ifelse(diagnosis=="HD","HD","MS")) %>% 
  mutate(sex = ifelse(sex %in% c("F","Female"),"Female","Male"))

soma_untreated$diff <- soma_untreated$csf_protein_total - soma_untreated$csf_albumin
soma_untreated$csf_protein_total <- z_score(soma_untreated$csf_protein_total,log=TRUE)
soma_untreated$csf_albumin <- z_score(soma_untreated$csf_albumin,log=TRUE)
soma_untreated$diff <- z_score(soma_untreated$diff,log=TRUE)

soma_untreated_adj <- left_join(soma_untreated_adj,model_predictions,by=c("sampleid","patientcode"))

############################################################
############## Examining accelerated aging #################
############################################################
soma_untreated <- soma_untreated %>%

  mutate(combined_diagnosis = ifelse(diagnosis=="HD","HD",ifelse(diagnosis == "RR-MS","RR-MS","P-MS"))) %>%
  mutate(combined_diagnosis = factor(combined_diagnosis,levels = c("HD","RR-MS","P-MS"))) %>%
  # mutate(diagnosis = ifelse(diagnosis == "RR-MS",))
  mutate(diagnosis = factor(diagnosis,levels = c("HD","RR-MS","SP-MS","PP-MS"))) %>%
  mutate(ms_diagnosis = factor(ms_diagnosis,levels=c("HD","MS")))

# soma_untreated %>% 
#   ggplot(aes(x=diagnosis,y=age)) +
#   geom_boxplot()
# 
# age_detected <- read_lines("./results/age_sex_data/age_selected_hd_detected.txt")
# age_formula <- as.formula(paste("age ~",paste(age_detected,collapse = "+")))
# 
library(glmnet)
library(broom)

hd_mean <- sapply(soma_untreated[soma_untreated$diagnosis=="HD",new_markers],mean)
hd_sd <- sapply(soma_untreated[soma_untreated$diagnosis=="HD",new_markers],sd)

x <- scale(as.matrix(soma_untreated[soma_untreated$diagnosis=="HD",new_markers]))
y <- soma_untreated$age[soma_untreated$diagnosis=="HD"]
y_mean <- mean(y)
y <- y - mean(y)

x_test <- scale(as.matrix(soma_untreated[soma_untreated$diagnosis!="HD",new_markers]),
                center=hd_mean,scale = hd_sd)
y_test <- soma_untreated$age[soma_untreated$diagnosis!="HD"]
y_test <- y_test - y_mean

set.seed(1234)
K <- 5
foldid <- sample(rep(1:K,ceiling(dim(x)[1]/K)),dim(x)[1],replace=FALSE)
alpha_seq <- seq(0.1,1,length=10)
for(i in alpha_seq){
  # set.seed(1234)
  cv_out <- cv.glmnet(x,y,alpha=i,standardize=FALSE,intercept=FALSE,foldid = foldid)
  print(cv_out$cvm[which(cv_out$lambda==cv_out$lambda.min)])
  # print(cv_out$cvm[which(cv_out$lambda==cv_out$lambda.1se)])
}
# 0.7 and 0.9
set.seed(1234)
cv_out <- cv.glmnet(x,y,alpha=0.9,standardize=FALSE,intercept=FALSE,nfolds = 5)
plot(cv_out)
plot(cv_out$glmnet.fit,"lambda")

glmnet_out <- tidy(cv_out$glmnet.fit) 
glmnet_fit_best <- filter(glmnet_out,lambda==cv_out$lambda.1se)
# glmnet_fit_best <- filter(glmnet_out,lambda==cv_out$lambda.min)
glmnet_fit_best <- glmnet_fit_best %>%
  select(term,estimate) %>%
  arrange(desc(abs(estimate)))
glmnet_fit_best <- as_tibble(cbind(clean_marker(glmnet_fit_best[[1]]),glmnet_fit_best[,-1]))
# write_csv(glmnet_fit_best,"elastic_net.csv")

# cor.test(y_test,y_test_pred <- predict(cv_out$glmnet.fit,x_test,s=cv_out$lambda.1se))

# hist(y_test-y_test_pred)
# boxplot(as.numeric(y_test-y_test_pred) ~ soma_untreated$diagnosis[soma_untreated$diagnosis!="HD"])


soma_untreated$predicted_age <- predict(cv_out$glmnet.fit,scale(as.matrix(soma_untreated[,new_markers]),center=hd_mean,scale = hd_sd),s=cv_out$lambda.1se) + y_mean
# soma_untreated$predicted_age <- predict(cv_out$glmnet.fit,scale(as.matrix(soma_untreated[,new_markers]),center=hd_mean,scale = hd_sd),s=cv_out$lambda.min) + y_mean
soma_untreated$residual_age <- with(soma_untreated,age-predicted_age)

soma_plot <- soma_untreated %>% 
  mutate(diagnosis = factor(diagnosis,levels=c("HD","RR-MS","SP-MS","PP-MS"))) %>% 
  mutate(ms_diagnosis = as.factor(ms_diagnosis))

p1 <- soma_plot %>% 
  rename(Diagnosis = diagnosis) %>% 
  ggplot(aes(x=predicted_age,y=age,group=Diagnosis,color=Diagnosis)) +
  geom_point() +
  geom_smooth(method="lm") +
  # geom_abline(intercept=0,slope=1,size=1.5) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

p2 <- soma_plot %>% 
  mutate(diagnosis = as.character(diagnosis)) %>% 
  mutate(diagnosis = str_replace_all(diagnosis,"-","")) %>% 
  mutate(diagnosis = factor(ifelse(diagnosis=="HD","HV",diagnosis),levels=c("HV","RRMS","SPMS","PPMS"))) %>% 
  rename(Diagnosis = diagnosis) %>% 
  mutate(residual_age = -1*residual_age) %>% 
  ggplot(aes(x=Diagnosis,y=residual_age,color=Diagnosis)) +
  geom_boxplot() +
  # geom_jitter(height=0,width=0.1) +
  geom_beeswarm() +
  # geom_text(aes(x=1,y=20,label="Overpredicting age"),color="black",size=4) +
  # geom_text(aes(x=1,y=-20,label="Underpredicting age"),color="black",size=4) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        legend.position="none",
        panel.grid = element_blank()) 

# tiff(my_filename("./results/elasticnet/elasticnet_slopes.tiff"),units="in",res=300,width=5,height=4.35)
# grid.arrange(p1)
# dev.off()

tiff(my_filename("./results/elasticnet/elasticnet_residuals.tiff"),units="in",res=300,width=5.8,height=4.8)
grid.arrange(p2)
dev.off()

# tiff("examining_excelarated_aging.tiff",units="in",res=300,width=14,height=5)
# grid.arrange(p1,p3,ncol=2)
# dev.off()
# 
# par(mfrow=c(1,2))
# boxplot(age~diagnosis,data=soma_plot,ylim=range(with(soma_plot,c(age,predicted_age))),
#         ylab="Observed Age",main="Observed Ages")
# boxplot(predicted_age~diagnosis,data=soma_plot,ylim=range(with(soma_plot,c(age,predicted_age))),
#         ylab="Predicted Age",main="Predicted Ages")
# par(mfrow=c(1,1))


soma_plot %>% 
  gather(key,value,age,predicted_age) %>% 
  ggplot(aes(x=key,y=value,group=patientcode)) +
  facet_wrap(~diagnosis) +
  geom_point() +
  geom_line()

Anova(lm(age~predicted_age,data=soma_plot,subset=diagnosis=="HD"))
summary(lm(age~predicted_age,data=soma_plot,subset=diagnosis=="HD"))

with(subset(soma_plot,diagnosis=="HD"),sqrt(mean((age - predicted_age)^2)))

Anova(lm(age~predicted_age*diagnosis,data=soma_plot))
summary(lm(age~predicted_age*diagnosis,data=soma_plot))

Anova(lm(age~predicted_age*ms_diagnosis,data=soma_plot))
summary(lm(age~predicted_age*ms_diagnosis,data=soma_plot))

Anova(lm(residual_age~diagnosis,data=soma_plot))
summary(lm(residual_age~diagnosis,data=soma_plot))
cld(glht(lm(residual_age~diagnosis,data=soma_plot),linfct=mcp(diagnosis="Tukey")))

Anova(lm(residual_age~ms_diagnosis,data=soma_plot))
summary(lm(residual_age~ms_diagnosis,data=soma_plot))
cld(glht(lm(residual_age~ms_diagnosis,data=soma_plot),linfct=mcp(ms_diagnosis="Tukey")))
# 
# str(cv_out$glmnet.fit$beta[,which])
# 
age_detected <- read_lines("./results/age_sex_data/age_selected_hd_detected.txt")
age_formula <- as.formula(paste("age ~",paste(age_detected,collapse = "+")))
# age_formula <- as.formula(paste("age ~",paste(age_to_adjust,collapse = "+")))


Anova(hd_mod <- lm(age_formula,data=soma_untreated,subset = diagnosis == "HD"))
Anova(hd_mod <- update(hd_mod,.~.-SL011628))
Anova(hd_mod <- update(hd_mod,.~.-SL007631))
Anova(hd_mod <- update(hd_mod,.~.-SL005233))
Anova(hd_mod <- update(hd_mod,.~.-SL001995))
Anova(hd_mod <- update(hd_mod,.~.-SL007804))
Anova(hd_mod <- update(hd_mod,.~.-SL010612))
Anova(hd_mod <- update(hd_mod,.~.-SL012395))
Anova(hd_mod <- update(hd_mod,.~.-SL010388)) # Stop here if 0.05 cutoff level
# Anova(hd_mod <- update(hd_mod,.~.-SL006230))
# 
# Anova(hd_mod <- update(hd_mod,.~.-SL004811))
# Anova(hd_mod <- update(hd_mod,.~.-SL000047)) # Stop here if 0.01 cutoff level
# 
# summary(hd_mod <- lm(age~SL009324,data=soma_untreated,subset = diagnosis == "HD"))
# summary(hd_mod <- lm(age~SL009324+SL003869,data=soma_untreated,subset = diagnosis == "HD"))
# summary(hd_mod <- lm(age~SL009324+SL003323+SL003869,data=soma_untreated,subset = diagnosis == "HD"))
# 
# summary(hd_mod <- lm(age~SL003323+SL003869,data=soma_untreated,subset = diagnosis == "HD"))
# summary(hd_mod <- lm(age~SL003323+SL003869+SL009324,data=soma_untreated,subset = diagnosis == "HD"))
# 
# library(psych)
# pairs.panels(subset(soma_untreated[,c("age","SL003323","SL003869","SL009324")]))
# 
# plot(hd_mod,which=1)
# 
# 
soma_untreated$predicted_age <- predict(hd_mod,newdata=soma_untreated)
soma_untreated$age_resid <- with(soma_untreated,age - predicted_age)
lims <- with(soma_untreated,range(age,predicted_age))
soma_untreated %>%
  ggplot(aes(x=predicted_age,y=age)) +
  facet_wrap( ~diagnosis) +
  geom_point() +
  geom_abline(slope=1,intercept = 0) +
  xlab("Predicted Age") +
  ylab("Observed Age") +
  theme_bw() +
  ylim(lims) +
  xlim(lims)

soma_untreated %>%
  ggplot(aes(x=diagnosis,y=age_resid)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(height=0,width=0.1) +
  xlab("") +
  ylab("Age - Predicted Age") +
  ggtitle("Excelerated Aging Residuals") +
  theme_bw()
Anova(full_mod <- lm(age_resid~diagnosis,data=soma_untreated))
cld(glht(full_mod,linfct=mcp(diagnosis="Tukey")))
# Anova(full_mod <- lm(age_resid~combined_diagnosis,data=soma_untreated))
# cld(glht(full_mod,linfct=mcp(combined_diagnosis="Tukey")))
# 
# plot(full_mod,which=1)
# plot(allEffects(full_mod))
# 
# detected_changes <- sort(sapply(age_detected,function(soma_char){
#   # soma_char <- "SL003869"
#   mod <- lm(soma_untreated[[soma_char]]~age*ms_diagnosis,data=soma_untreated)
#   return(summary(mod)$coefficients[4,4])
# }))
# # detected_elevations <- sort(sapply(age_detected,function(soma_char){
# #   # soma_char <- "SL003869"
# #   mod <- lm(soma_untreated[[soma_char]]~age+ms_diagnosis,data=soma_untreated)
# #   return(summary(mod)$coefficients[3,4])
# # }))
# 
# soma_untreated %>% 
#   ggplot(aes(x=age,y=SL003869,group=ms_diagnosis,color=ms_diagnosis)) +
#   geom_point() +
#   geom_smooth(method="lm") +
#   xlab("Age") +
#   ylab("GDF15")
# 
# soma_untreated %>% 
#   ggplot(aes(x=age,y=SL002785,group=ms_diagnosis,color=ms_diagnosis)) +
#   geom_point() +
#   geom_smooth(method="lm") +
#   xlab("Age") +
#   ylab("NPPB")
# 
# soma_untreated %>% 
#   ggplot(aes(x=age,y=SL003323,group=ms_diagnosis,color=ms_diagnosis)) +
#   geom_point() +
#   geom_smooth(method="lm") +
#   xlab("Age") +
#   ylab("CCL18")
# 
############################################################
############## Loading CNS Expression Dat ##################
############################################################
expression_dat <- read_excel("./data/1-s2.0-S0896627315010193-mmc3.xlsx",sheet="Human data only",
                             skip=3,col_names = FALSE)

expression_dat2 <- read_excel("./data/1-s2.0-S0896627315010193-mmc3.xlsx",sheet="Human and mouse",
                             skip=3,col_names = FALSE)[,1:42]

names(expression_dat) <- c("gene",paste("V",1:41,sep=""))
names(expression_dat2) <- c("gene",paste("V",1:41,sep=""))

expression_dat <- bind_rows(expression_dat,expression_dat2) %>% 
  gather(key,expression,-gene) %>% 
  mutate(gene = toupper(gene)) %>% 
  mutate(gene = ifelse(gene == "C10ORF54","C10orf54",gene))
remove(expression_dat2)

nice_names <- c(rep("GBM/peri-tumor astrocytes",4),
                rep("Sclerotic hippocampi astrocytes",4),
                rep("Fetal astrocytes",6),
                rep("Mature astrocytes",12),
                "Neurons",
                rep("Oligodendrocytes",5),
                rep("Microglia/Macrophage",3),
                rep("Endothelial",2),
                rep("Whole cortex",4))
name_dat <- tibble(key = paste("V",1:41,sep=""),
                   cell_type = nice_names)

# levs <- c("Mature astrocytes","Neurons","Oligodendrocytes","Microglia/Macrophage","Endothelial")
# levs <- c("GBM/peri-tumor astrocytes","Sclerotic hippocampi astrocytes","Fetal astrocytes",
#           "Mature astrocytes","Neurons","Oligodendrocytes","Microglia/Macrophage","Endothelial",
#           "Whole cortex")
# levs <- c("GBM/peri-tumor astrocytes","Sclerotic hippocampi astrocytes",
#           "Mature astrocytes","Neurons","Oligodendrocytes","Microglia/Macrophage","Endothelial",
#           "Whole cortex")

levs <- c("Mature astrocytes","Neurons","Oligodendrocytes","Microglia/Macrophage","Endothelial",
          "Whole cortex")


expression_dat <- as_tibble(merge(expression_dat,name_dat,by="key"))
expression_dat <- expression_dat %>% 
  filter(cell_type %in% levs) %>% 
  mutate(cell_type = factor(cell_type,levels=levs)) %>% 
  select(-key) %>% 
  group_by(gene,cell_type) %>% 
  summarize_all(funs(mean,sd,length)) %>% 
  mutate(se = ifelse(is.nan(sd),0,sd/sqrt(length))) %>% 
  select(-sd,-length)

tomatch <- read_excel("./data/translation file.xlsx",col_names = FALSE)
tomatch <- tomatch[,2:4]
names(tomatch) <- c("soma","gene","protein")

fixed_gene <- read_excel("./data/missing_gene.xlsx")

tomatch <- merge(tomatch,fixed_gene,by="gene",all.x=TRUE) %>% 
  mutate(gene = ifelse(is.na(new_name),gene,new_name)) %>% 
  select(-new_name) %>% 
  as_tibble()

# tomatch_gene <- tomatch %>% 
#   select(gene) %>% 
#   remove_dups()
# 
# gene_match <- left_join(tomatch_gene,expression_dat,by="gene")
# gene_missing <- gene_match %>% 
#   filter(is.na(cell_type)) %>% 
#   select(gene) %>% 
#   remove_dups()
# 
# write_csv(gene_missing,"./data/missing_gene.csv")

expression_dat <- dplyr::inner_join(tomatch,expression_dat,by="gene")
expression_dat_raw <- expression_dat %>%
  mutate(mean = log(mean))

expression_dat <- expression_dat %>%
  mutate(mean = ifelse(mean == 0.1,NA,mean)) %>%
  group_by(cell_type) %>%
  mutate(mean = z_score(mean,log=TRUE)) %>%
  mutate(mean = ifelse(is.na(mean),min(mean,na.rm=TRUE),mean)) %>%
  ungroup()

expression_dat <- expression_dat %>% 
  select(-se) %>% 
  spread(cell_type,mean) 

##########################################################################
################## Pub Dat ###############################################
##########################################################################
pub_dat <- read_excel("./data/Genomic atlas of serum proteins Tables.xlsx",sheet=2,skip=4,
                      col_names = FALSE,na=c("",NA,"--"))

names(pub_dat) <- c("id","protein","target","uniprot","r2_adj","sub_pass","age_beta_pub","age_se_pub","age_pval_pub",
                    "female_beta_pub","female_se_pub","female_pval_pub","bmi_beta_pub","bmi_se_pub","bmi_pval_pub",
                    "egfr_beta_pub","egfr_se_pub","egfr_pval_pub")

pub_dat <- pub_dat[,-c(19)]

pub_dat <- pub_dat %>% 
  mutate(age_pval_pub = 10^(-1*age_pval_pub),
         female_pval_pub = 10^(-1*female_pval_pub)) %>% 
  filter(age_pval_pub <= 1e-5 | female_pval_pub<=1e-5) %>%
  select(id:female_pval_pub) %>% 
  select(-protein,-target,-uniprot,-r2_adj,-sub_pass)

tomatch <- read_excel("translation1.3k.xlsx")
tomatch <- tomatch %>% 
  rename(gene = EntrezGeneSymbol,protein=Target,protein_full=TargetFullName,marker=SomaId,id=SeqId) %>% 
  mutate(id = paste(gene,".",id,sep="")) %>% 
  mutate(id = str_replace_all(id,"-",".")) %>% 
  mutate(id = str_replace_all(id,"_",".")) %>% 
  mutate(id = str_replace_all(id,",",".")) %>%
  mutate(id = str_replace_all(id,"@",".")) %>%
  mutate(id = str_replace_all(id," ",".")) %>%
  mutate(id = str_replace_all(id,"(\\.)(\\.)",".")) %>%
  select(-UniProt,-EntrezGeneID)

tomatch <- left_join(tomatch,pub_dat,by="id")
tomatch <- tomatch %>% 
  rename(soma = marker)

# overlap <- function(...){
#   ll <- list(...)
#   n <- length(ll)
#   ll <- table(unlist(ll))
#   ll <- names(ll)[which(ll==n)]
#   return(ll)
# }


# set.seed(1324)
# x <- rnorm(100)
# y <- 2 + 1 * x + rnorm(100,0,1)
# 
# plot(y~x)
# abline(mod <- lm(y~x))
# plot(I(y+predict(mod <- lm(y~x))) ~ x)
# 
# summary(mod)
# summary(lm((I(y+predict(mod)) ~ x)))

##########################################################################
################## Sex - Unadjusted ######################################
##########################################################################
sex <- read_excel("./data/Data for String and Zhang.xlsx",sheet = "sex")
sex_disc <- merge(sex,expression_dat,by="soma")
sex_dir <- sex_disc %>% 
  filter(soma %in% sex_disc$soma) 

sex_disc <- sex_disc %>% 
  filter(soma %in% sex_dir$soma) %>% 
  left_join(select(tomatch,-protein_full,-protein,-gene),by=c("soma"))

sex_dir <- sex_dir %>% 
  .$direction

sex_neuro <- sex_disc[,levs] %>% 
  as.matrix()

sex_soma <- t(soma_untreated[,sex_disc$soma])
row.names(sex_soma) <- sex_disc$gene
row.names(sex_neuro) <- sex_disc$gene

draw(Heatmap(sex_neuro,clustering_method_rows = "ward.D2",cluster_columns = TRUE,
             clustering_distance_rows = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = FALSE,
             show_row_names = TRUE,column_title = "CNS Expression",width = 1,
             col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),row_names_side = "left",
             # col=colorRamp2(c(0.1,1,10),c("green","black","red"),space="sRGB"),row_names_side = "left",
             split = sex_dir,bottom_annotation = HeatmapAnnotation(
               text = anno_text(levs,rot=45,offset = unit(1, "npc"), just = "right"),
               annotation_height = max_text_width(levs)
             ),
             heatmap_legend_param = list(legend_direction = "horizontal",title_position="leftcenter")) +
       Heatmap(sex_soma[,soma_untreated$sex=="Male"],clustering_method_columns = "ward.D2",cluster_columns = TRUE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "Male",width=4.5,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               heatmap_legend_param = list(legend_direction = "horizontal",title_position="leftcenter")) + 
       Heatmap(sex_soma[,soma_untreated$sex=="Female"],clustering_method_columns = "ward.D2",cluster_columns = TRUE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "Female",width=5.5,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               heatmap_legend_param = list(legend_direction = "horizontal",title_position="leftcenter")),
     column_title = "Proteins with Sex Associations",
     show_heatmap_legend=TRUE,heatmap_legend_side="bottom")

sex_dend <- as.dendrogram(hclust(as.dist(1-cor(sex_soma)),method="ward.D2"))


rowAnnotation(data.frame(Beta = sex_disc$female_beta_pub),
              col = list(Beta = colorRamp2(c(-0.04,0,.04), c("green","white", "red"),space="sRGB")))

# tiff(my_filename("./results/age_sex_adjust_fig/Gender_Heatmap.tiff"),units="in",res=300,width=6.666,height=7.5)
tiff(my_filename("./results/age_sex_adjust_fig/Gender_Heatmap.tiff"),units="in",res=300,width=11.25,height=7.5)
draw(Heatmap(sex_soma[,soma_untreated$sex=="Male"],clustering_method_columns = "ward.D2",
             cluster_columns = TRUE,cluster_rows=TRUE,
             clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
             show_row_names = TRUE,column_title = "Males",width=4.5,
             split = sex_dir,row_names_side = "left",
             col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
             heatmap_legend_param = list(legend_direction = "horizontal",title_position="leftcenter")) + 
       Heatmap(sex_soma[,soma_untreated$sex=="Female"],clustering_method_columns = "ward.D2",
               cluster_columns = TRUE,cluster_rows=TRUE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "Females",width=5.5,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               heatmap_legend_param = list(legend_direction = "horizontal")),
     column_title = "Proteins with Sex Associations",
     show_heatmap_legend=TRUE,heatmap_legend_side="bottom",gap=unit(1,"mm"))
dev.off()

tiff(my_filename("./results/age_sex_adjust_fig/Gender_Heatmap.tiff"),units="in",res=300,width=6.666,height=7.5)
draw(rowAnnotation("INTERVAL Female Beta" = sex_disc$female_beta_pub,
                   col = list("INTERVAL Female Beta" = colorRamp2(c(-1.4,0,1.4), c("green","white", "red"),space="sRGB"))) + 
       Heatmap(sex_soma[,soma_untreated$sex=="Male"],clustering_method_columns = "ward.D2",
             cluster_columns = TRUE,cluster_rows=TRUE,
             clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
             show_row_names = TRUE,column_title = "Males",width=4.5,
             split = sex_dir,row_names_side = "left",
             col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
             heatmap_legend_param = list(legend_direction = "horizontal",title_position="leftcenter")) + 
       Heatmap(sex_soma[,soma_untreated$sex=="Female"],clustering_method_columns = "ward.D2",
               cluster_columns = TRUE,cluster_rows=TRUE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "Females",width=5.5,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               heatmap_legend_param = list(legend_direction = "horizontal")),
     column_title = "Proteins with Sex Associations",row_dend_side = "left",
     show_heatmap_legend=TRUE,heatmap_legend_side="bottom",gap=unit(1,"mm"),row_sub_title_side = "left")
dev.off()


######################################################################
################## Age Associations ##################################
######################################################################
# age <- read_excel("./data/Data for String and Zhang.xlsx",sheet = "age")
# 
# ####################
# totprot_cor <- tibble(soma = age$soma) %>% 
#   mutate(total_protein_cor = sapply(soma,function(x){cor.test(soma_untreated[[x]],soma_untreated[["csf_protein_total"]],method="spearman")$est})) %>% 
#   mutate(total_protein_pval = sapply(soma,function(x){cor.test(soma_untreated[[x]],soma_untreated[["csf_protein_total"]],method="spearman")$p.value})) %>% 
#   mutate(albumin_cor = sapply(soma,function(x){cor.test(soma_untreated[[x]],soma_untreated[["csf_albumin"]],method="spearman")$est})) %>% 
#   mutate(albumin_pval = sapply(soma,function(x){cor.test(soma_untreated[[x]],soma_untreated[["csf_albumin"]],method="spearman")$p.value})) %>% 
#   mutate(diff_cor = sapply(soma,function(x){cor.test(soma_untreated[[x]],soma_untreated[["csf_protein_total"]]-soma_untreated[["csf_albumin"]],method="spearman")$est})) %>% 
#   mutate(diff_pval = sapply(soma,function(x){cor.test(soma_untreated[[x]],soma_untreated[["csf_protein_total"]]-soma_untreated[["csf_albumin"]],method="spearman")$p.value})) %>% 
#   arrange(diff_pval)
# 
# soma_untreated_ms <- filter(soma_untreated,diagnosis !="HD")
# soma_untreated_hd <- filter(soma_untreated,diagnosis =="HD")
# 
# totprot_cor_ms <- tibble(soma = str_subset(names(soma_untreated_ms),"^SL")) %>% 
#   mutate(age_adjusted = ifelse(soma %in% age$soma,"Yes","No")) %>% 
#   mutate(pipeline_marker = ifelse(soma %in% unique(unlist(strsplit(pipeline_rats,"/"))),"Yes","No")) %>% 
#   mutate(total_protein_cor = sapply(soma,function(x){cor.test(soma_untreated_ms[[x]],soma_untreated_ms[["csf_protein_total"]],method="spearman")$est})) %>% 
#   mutate(total_protein_pval = sapply(soma,function(x){cor.test(soma_untreated_ms[[x]],soma_untreated_ms[["csf_protein_total"]],method="spearman")$p.value})) %>% 
#   mutate(albumin_cor = sapply(soma,function(x){cor.test(soma_untreated_ms[[x]],soma_untreated_ms[["csf_albumin"]],method="spearman")$est})) %>% 
#   mutate(albumin_pval = sapply(soma,function(x){cor.test(soma_untreated_ms[[x]],soma_untreated_ms[["csf_albumin"]],method="spearman")$p.value})) %>% 
#   mutate(diff_cor = sapply(soma,function(x){cor.test(soma_untreated_ms[[x]],soma_untreated_ms[["csf_protein_total"]]-soma_untreated_ms[["csf_albumin"]],method="spearman")$est})) %>% 
#   mutate(diff_pval = sapply(soma,function(x){cor.test(soma_untreated_ms[[x]],soma_untreated_ms[["csf_protein_total"]]-soma_untreated_ms[["csf_albumin"]],method="spearman")$p.value})) %>% 
#   arrange(diff_pval)
# 
# totprot_cor_hd <- tibble(soma = str_subset(names(soma_untreated_hd),"^SL")) %>%
#   mutate(age_adjusted = ifelse(soma %in% age$soma,"Yes","No")) %>%
#   mutate(pipeline_marker = ifelse(soma %in% unique(unlist(strsplit(pipeline_rats,"/"))),"Yes","No")) %>%
#   mutate(total_protein_cor = sapply(soma,function(x){cor.test(soma_untreated_hd[[x]],soma_untreated_hd[["csf_protein_total"]],method="spearman")$est})) %>%
#   mutate(total_protein_pval = sapply(soma,function(x){cor.test(soma_untreated_hd[[x]],soma_untreated_hd[["csf_protein_total"]],method="spearman")$p.value})) %>%
#   mutate(albumin_cor = sapply(soma,function(x){cor.test(soma_untreated_hd[[x]],soma_untreated_hd[["csf_albumin"]],method="spearman")$est})) %>%
#   mutate(albumin_pval = sapply(soma,function(x){cor.test(soma_untreated_hd[[x]],soma_untreated_hd[["csf_albumin"]],method="spearman")$p.value})) %>%
#   mutate(diff_cor = sapply(soma,function(x){cor.test(soma_untreated_hd[[x]],soma_untreated_hd[["csf_protein_total"]]-soma_untreated_hd[["csf_albumin"]],method="spearman")$est})) %>%
#   mutate(diff_pval = sapply(soma,function(x){cor.test(soma_untreated_hd[[x]],soma_untreated_hd[["csf_protein_total"]]-soma_untreated_hd[["csf_albumin"]],method="spearman")$p.value})) %>%
#   arrange(diff_pval)
# 
# # clean_marker(age$soma) %>% 
# #   rename(soma = marker) %>% 
# #   left_join(totprot_cor_hd) %>% 
# #   select(-age_adjusted,-pipeline_marker) %>% 
# #   arrange(total_protein_pval) %>% 
# #   write_csv("./results/Age adjusted somamer correlations with CSF protein - HD - 20190419.csv")
# # clean_marker(age$soma) %>% 
# #   rename(soma = marker) %>% 
# #   left_join(totprot_cor_ms) %>% 
# #   select(-age_adjusted,-pipeline_marker) %>% 
# #   arrange(total_protein_pval) %>% 
# #   write_csv("./results/Age adjusted somamer correlations with CSF protein - MS - 20190419.csv")
# 
# # totprot_cor_hd <- arrange(totprot_cor_hd,soma)
# # totprot_cor_ms <- arrange(totprot_cor_ms,soma)
# # plot(totprot_cor_hd$total_protein_cor[totprot_cor_hd$age_adjusted=="Yes"],
# #      totprot_cor_ms$total_protein_cor[totprot_cor_ms$age_adjusted=="Yes"])
# # abline(0,1)
# # grid.arrange(totprot_cor_ms %>% 
# #                gather(key,value,total_protein_cor,albumin_cor,diff_cor) %>% 
# #                mutate(key = ifelse(key =="diff_cor","Total - Albumin",ifelse(key=="albumin_cor","CSF Albumin","CSF Protein"))) %>% 
# #                mutate(key = factor(key,levels=c("CSF Protein","CSF Albumin", "Total - Albumin"))) %>% 
# #                ggplot(aes(x=age_adjusted,y=value)) +
# #                facet_wrap(~key,scales="fixed") +
# #                geom_boxplot() +
# #                geom_jitter(alpha=0.4,height=0,width=0.1) +
# #                xlab("Age adjusted?") +
# #                ylab("Spearman Correlation") +
# #                ggtitle("Comparing all somamers and those that were age adjusted") +
# #                theme_bw()+
# #                theme(plot.title = element_text(hjust=0.5)),
# #              totprot_cor_ms %>% 
# #                gather(key,value,total_protein_cor,albumin_cor,diff_cor) %>% 
# #                mutate(key = ifelse(key =="diff_cor","Total - Albumin",ifelse(key=="albumin_cor","CSF Albumin","CSF Protein"))) %>% 
# #                mutate(key = factor(key,levels=c("CSF Protein","CSF Albumin", "Total - Albumin"))) %>% 
# #                ggplot(aes(x=pipeline_marker,y=value)) +
# #                facet_wrap(~key,scales="fixed") +
# #                geom_boxplot() +
# #                geom_jitter(alpha=0.4,height=0,width=0.1) +
# #                xlab("Selected in pipeline?") +
# #                ylab("Spearman Correlation")+
# #                ggtitle("Comparing all somamers and those that pipeline selected") +
# #                theme_bw()+
# #                theme(plot.title = element_text(hjust=0.5)),ncol=1)
# 
# 
# # pairs.panels(soma_untreated[,c("csf_protein_total","csf_albumin","diff")],
# #              ellipse = FALSE,smooth = FALSE,lm=TRUE,stars=FALSE)
# 
####################
age <- read_excel("./data/Data for String and Zhang.xlsx",sheet = "age")
age_disc <- merge(age,expression_dat,by="soma")

age_dir <- age_disc %>%
  filter(soma %in% age_disc$soma) %>%
  mutate(dir = paste(direction,direction2,sep="-"))

age_disc <- age_disc %>% 
  filter(soma %in% age_dir$soma) %>% 
  left_join(select(tomatch,-protein_full,-protein,-gene),by=c("soma"))

age_neuro <- age_disc[which(age_disc$soma %in% age_dir$soma),levs] %>% 
  as.matrix()

soma_untreated <- soma_untreated %>% 
  arrange(age)
age_soma <- t(soma_untreated[,age_dir$soma])

row.names(age_soma) <- age_disc$gene[which(age_disc$soma %in% age_dir$soma)]
row.names(age_neuro) <- age_disc$gene[which(age_disc$soma %in% age_dir$soma)]
age_detected <- with(age_dir,soma[direction2 == "Detected"])
row_cols <- ifelse(age_dir$soma %in% age_detected,"maroon","black")
age_dir <- age_dir %>%
  .$direction

# age_cordat <- age_soma
# age_cordat <- rbind(soma_untreated$age,age_soma)[c(1,1+row_order(age_severity)[[1]],1+row_order(age_severity)[[2]]),]
# row.names(age_cordat)[1] <- "Age"
# cor_col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
#                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
#                                   "#4393C3", "#2166AC", "#053061")))(200)
# # cor_col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
# #                                   "#92C5DE",
# #                                   "#4393C3", "#2166AC", "#053061")))(200)
# corrplot(cor(t(age_cordat[,soma_untreated$diagnosis!="HD"])),
#          method="square",tl.pos = "l",addgrid.col = NA,type="lower",
#          tl.cex=0.45,title="",
#          mar = c(0, 0, 2, 0),col=cor_col,cl.pos="b")
# corrplot(cor(t(age_cordat[,soma_untreated$diagnosis=="HD"])),
#          method="square",tl.pos = "t",addgrid.col = NA,type="upper",
#          tl.cex=0.45,
#          mar = c(0, 0, 2, 0),col=cor_col,add=TRUE,cl.pos="n")
# 

row_order(age_severity)[[1]]

draw(age_severity <- rowAnnotation(data.frame(Beta = age_disc$age_beta_pub),
                   col = list(Beta = colorRamp2(c(-0.04,0,.04), c("green","white", "red"),space="sRGB"))) + 
       Heatmap(age_neuro,clustering_method_rows = "ward.D2",cluster_columns = TRUE,
               clustering_distance_rows = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = FALSE,
               show_row_names = TRUE,column_title = "CNS Expression",width = 1,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               row_names_side = "left",
               split = age_dir, 
               # combined_name_fun = NULL,
               row_names_gp = gpar(col = row_cols, fontsize = 7),show_column_dend = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(levs,rot=45,offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(levs)
               )) +
       Heatmap(age_soma[,soma_untreated$diagnosis=="HD"],clustering_method_rows = "ward.D2",cluster_columns = FALSE,
               clustering_distance_rows = "spearman",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "HD Patient Expression",width=1.8,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               row_names_gp = gpar(col = row_cols, fontsize = 7),
               row_names_side = "left",split = age_dir,
               top_annotation = HeatmapAnnotation(
                 # data.frame(age = soma_untreated[soma_untreated$diagnosis=="HD",][["age"]]),
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis=="HD",][["age"]],axis=TRUE,side="left"),
                 "Total Protein" =soma_untreated[soma_untreated$diagnosis=="HD",][["csf_protein_total"]],
                 "Total - Albumin" = soma_untreated[soma_untreated$diagnosis=="HD",][["diff"]],
                 col = list("Total Protein" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
                            "Total - Albumin" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
                 show_annotation_name = TRUE,
                 annotation_name_offset = unit(1, "cm"),
                 annotation_name_rot = c(0, 0, 0),
                 height = unit(2.25, "cm"),
                 annotation_name_side = "left",
                 annotation_height = c(3,0.75,0.75)),
               heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(5, "cm"), title_position = "lefttop")) +
       Heatmap(age_soma[,soma_untreated$diagnosis!="HD"],clustering_method_columns = "ward.D2",cluster_columns = FALSE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = TRUE,column_title = "MS Patient Expression",width=8.2,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               row_names_gp = gpar(col = row_cols, fontsize = 7),
               row_names_side = "right",split = age_dir,
               top_annotation = HeatmapAnnotation(
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis!="HD",][["age"]],
                                   axis=FALSE,axis_side="right"),
                 "Total Protein" =soma_untreated[soma_untreated$diagnosis!="HD",][["csf_protein_total"]],
                 "Total - Albumin" = soma_untreated[soma_untreated$diagnosis!="HD",][["diff"]],
                 col = list("Total Protein" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
                            "Total - Albumin" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
                 show_annotation_name = FALSE,
                 annotation_name_offset = unit(1, "cm"),
                 annotation_name_rot = c(0, 0, 0),
                 height = unit(2.25, "cm"),
                 annotation_height = c(3,0.75,0.75))),
     column_title = "Proteins with Age Associations",
     show_heatmap_legend=TRUE,row_dend_side = "left",row_sub_title_side = "left",
     show_annotation_legend=FALSE,heatmap_legend_side="bottom",
     padding = unit(c(2, 2, 2, 30), "mm"),gap=unit(1,"mm"))

draw(rowAnnotation(data.frame(Beta = age_disc$age_beta_pub),
                   col = list(Beta = colorRamp2(c(-0.04,0,.04), c("green","white", "red"),space="sRGB"))) + 
       Heatmap(age_soma[,soma_untreated$diagnosis=="HD"],clustering_method_rows = "ward.D2",cluster_columns = FALSE,
               clustering_distance_rows = "spearman",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = TRUE,column_title = "HD Patient Expression",width=1.8,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               row_names_gp = gpar(col = row_cols, fontsize = 7),
               row_names_side = "left",split = age_dir,
               top_annotation = HeatmapAnnotation(
                 # data.frame(age = soma_untreated[soma_untreated$diagnosis=="HD",][["age"]]),
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis=="HD",][["age"]],axis=TRUE,side="left"),
                 "Total Protein" =soma_untreated[soma_untreated$diagnosis=="HD",][["csf_protein_total"]],
                 "Total - Albumin" = soma_untreated[soma_untreated$diagnosis=="HD",][["diff"]],
                 col = list("Total Protein" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
                            "Total - Albumin" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
                 show_annotation_name = TRUE,
                 annotation_name_offset = unit(1, "cm"),
                 annotation_name_rot = c(0, 0, 0),
                 height = unit(2.25, "cm"),
                 annotation_name_side = "left",
                 annotation_height = c(3,0.75,0.75)),
               heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(5, "cm"), title_position = "lefttop")) +
       Heatmap(age_soma[,soma_untreated$diagnosis!="HD"],clustering_method_columns = "ward.D2",cluster_columns = FALSE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "MS Patient Expression",width=8.2,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               top_annotation = HeatmapAnnotation(
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis!="HD",][["age"]],
                                   axis=FALSE,axis_side="right"),
                 "Total Protein" =soma_untreated[soma_untreated$diagnosis!="HD",][["csf_protein_total"]],
                 "Total - Albumin" = soma_untreated[soma_untreated$diagnosis!="HD",][["diff"]],
                 col = list("Total Protein" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
                            "Total - Albumin" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
                 show_annotation_name = FALSE,
                 annotation_name_offset = unit(1, "cm"),
                 annotation_name_rot = c(0, 0, 0),
                 height = unit(2.25, "cm"),
                 annotation_height = c(3,0.75,0.75)),
               heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(5, "cm"), title_position = "lefttop")) + 
       Heatmap(age_neuro,clustering_method_rows = "ward.D2",cluster_columns = TRUE,
               clustering_distance_rows = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = FALSE,
               show_row_names = TRUE,column_title = "CNS Expression",width = 1,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),row_names_side = "right",
               split = age_dir, 
               # combined_name_fun = NULL,
               row_names_gp = gpar(col = row_cols, fontsize = 7),show_column_dend = FALSE,
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(levs,rot=45,offset = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(levs)
               )),
     column_title = "Proteins with Age Associations",
     show_heatmap_legend=TRUE,row_dend_side = "left",row_sub_title_side = "left",
     show_annotation_legend=FALSE,heatmap_legend_side="bottom",
     padding = unit(c(2, 2, 2, 30), "mm"),gap=unit(1,"mm"))

hd_age_mat <- age_soma[,soma_untreated$diagnosis=="HD"]
ms_age_mat <- age_soma[,soma_untreated$diagnosis!="HD"]

# tiff(my_filename("./results/age_sex_adjust_fig/Age_Heatmap.tiff"),units="in",res=300,width=6.666,height=7.5)
# tiff(my_filename("./results/age_sex_adjust_fig/Age_Heatmap.tiff"),units="in",res=300,width=11.25,height=7.5)
tiff(my_filename("./results/age_sex_adjust_fig/Age_Heatmap.tiff"),units="in",res=300,width=11.75,height=7.5)
draw(age_severity <- Heatmap(hd_age_mat,clustering_method_rows = "ward.D2",cluster_columns = FALSE,
               clustering_distance_rows = "spearman",name="SomaScan Expression (Z-score)",
               show_column_names = TRUE,
               show_row_names = TRUE,column_title = "HD",width=1.8,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               row_names_gp = gpar(fontsize = 7),
               row_names_side = "left",split = age_dir,
               heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(5, "cm"), title_position = "lefttop"),
               top_annotation = HeatmapAnnotation(
                 # data.frame(age = soma_untreated[soma_untreated$diagnosis=="HD",][["age"]]),
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis=="HD",][["age"]],axis=TRUE,side="left"),
                 show_annotation_name = TRUE,
                 annotation_name_offset = unit(1, "cm"),
                 annotation_name_rot = c(0),
                 height = unit(1, "cm"),
                 annotation_name_side = "left",
                 annotation_height = c(3))) +
       Heatmap(ms_age_mat,clustering_method_columns = "ward.D2",cluster_columns = FALSE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "MS",width=8.2,
               row_names_gp = gpar(fontsize = 7), row_names_side = "right",
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"), na_col = "orange",
               top_annotation = HeatmapAnnotation(
                 # data.frame(age = soma_untreated[soma_untreated$diagnosis=="HD",][["age"]]),
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis!="HD",][["age"]],axis=FALSE,side="left"),
                 show_annotation_name = FALSE,
                 annotation_name_offset = unit(1, "cm"),
                 annotation_name_rot = c(0),
                 height = unit(1, "cm"),
                 annotation_name_side = "left",
                 annotation_height = c(3))),
     column_title = "Proteins with Age Associations",
     show_heatmap_legend=TRUE,row_dend_side = "left",row_sub_title_side = "left",
     heatmap_legend_side="bottom",annotation_legend_side="right",show_annotation_legend=FALSE,
     padding = unit(c(2, 2, 2, 30), "mm"),gap=unit(1,"mm"))
dev.off()

draw(age_severity <- rowAnnotation(data.frame(Beta = age_disc$age_beta_pub),
                                   col = list(Beta = colorRamp2(c(-0.04,0,.04), c("green","white", "red"),space="sRGB"))) + 
       Heatmap(hd_age_mat,clustering_method_rows = "ward.D2",cluster_columns = FALSE,
               clustering_distance_rows = "spearman",name="SomaScan Expression (Z-score)",
               show_column_names = TRUE,
               show_row_names = TRUE,column_title = "HD",width=1.8,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
               row_names_gp = gpar(col = row_cols, fontsize = 7), na_col = "orange",
               row_names_side = "left",split = age_dir,# row_gap = unit(5, "mm"),
               top_annotation = HeatmapAnnotation(
                 # data.frame(age = soma_untreated[soma_untreated$diagnosis=="HD",][["age"]]),
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis=="HD",][["age"]],axis=TRUE,side="left"),
                 "Total Protein" =soma_untreated[soma_untreated$diagnosis=="HD",][["csf_protein_total"]],
                 "Total - Albumin" = soma_untreated[soma_untreated$diagnosis=="HD",][["diff"]],
                 col = list("Total Protein" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
                            "Total - Albumin" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
                 show_annotation_name = TRUE,
                 annotation_name_offset = unit(1, "cm"),
                 annotation_name_rot = c(0, 0, 0),
                 height = unit(2.25, "cm"),
                 annotation_name_side = c("left","right","right"),
                 annotation_height = c(3,0.75,0.75),
                 show_legend = FALSE),
               heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(5, "cm"), title_position = "lefttop")) +
       Heatmap(ms_age_mat,clustering_method_columns = "ward.D2",cluster_columns = FALSE,
               clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
               show_row_names = TRUE,column_title = "MS",width=8.2,
               row_names_gp = gpar(col = row_cols, fontsize = 7), row_names_side = "right",
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"), na_col = "orange",
               top_annotation = HeatmapAnnotation(
                 Age = anno_points(soma_untreated[soma_untreated$diagnosis!="HD",][["age"]],
                                   axis=FALSE,axis_side="right"),
                 "Total Protein" =soma_untreated[soma_untreated$diagnosis!="HD",][["csf_protein_total"]],
                 "Total - Albumin" = soma_untreated[soma_untreated$diagnosis!="HD",][["diff"]],
                 col = list("Total Protein" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
                            "Total - Albumin" = colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
                 show_annotation_name = c(FALSE,TRUE,TRUE),
                 annotation_name_offset = unit(0.1, "cm"),
                 annotation_name_rot = c(0, 0, 0),na_col = "orange",
                 height = unit(2.25, "cm"),
                 annotation_name_side = c("left","right","right"),
                 annotation_height = c(3,0.75,0.75),
                 show_legend = FALSE)),
     column_title = "Proteins with Age Associations",
     show_heatmap_legend=TRUE,row_dend_side = "left",row_sub_title_side = "left",
     heatmap_legend_side="bottom",annotation_legend_side="right",show_annotation_legend=FALSE,
     padding = unit(c(2, 2, 2, 30), "mm"),gap=unit(1,"mm"))

# ,ht_gap = unit(c(5,5), "mm")

age_cordat <- age_soma
age_cordat <- rbind(soma_untreated$age,age_soma)[c(1,1+row_order(age_severity)[[1]],1+row_order(age_severity)[[2]]),]
row.names(age_cordat)[1] <- "Age"
cor_col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061")))(200)
# cor_col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
#                                   "#92C5DE",
#                                   "#4393C3", "#2166AC", "#053061")))(200)
corrplot(cor(t(age_cordat[,soma_untreated$diagnosis!="HD"])),
         method="square",tl.pos = "lt",addgrid.col = NA,type="lower",
         tl.cex=0.45,title="",
         tl.col=c("black",rep("red",53),rep("blue",22)),
         mar = c(0, 0, 2, 0),col=cor_col,cl.pos="b")
corrplot(cor(t(age_cordat[,soma_untreated$diagnosis=="HD"])),
         method="square",tl.pos = "n",addgrid.col = NA,type="upper",
         tl.cex=0.45,
         mar = c(0, 0, 2, 0),col=cor_col,add=TRUE,cl.pos="n")





##########################################################################
###################### Severity ##########################################
##########################################################################
# Clustering patients
pipeline_markers <- unique(unlist(strsplit(unlist(var_lists[c(8:10)]),"/")))
# pipeline_markers <- unique(unlist(strsplit(unlist(var_lists[c(3,4,7)]),"/")))
# pipeline_markers <- unique(unlist(strsplit(unlist(var_lists[c(3,5,7)]),"/")))
# pipeline_markers <- unique(unlist(strsplit(pipeline_rats,"/")))

severity_disc <- merge(tibble(soma=pipeline_markers),expression_dat,by="soma")
# severity_dir <- severity_disc %>% 
#   filter(soma %in% severity_disc$soma) %>% 
#   .$direction
severity_neuro <- severity_disc[,levs] %>% 
  as.matrix()

severity_data <- soma_untreated_adj %>% 
  filter(sampleid != "NDUS1187") %>%
  filter(diagnosis != "HD")

# col_mat <- severity_data[,pipeline_markers] %>% dist(method="manhattan")
col_mat <- severity_data[,pipeline_markers] %>% dist()

col_hc <-  col_mat %>% hclust(method="ward.D2")
# col_hc <-  col_mat %>% hclust(method="average")
col_dend <- col_hc %>% as.dendrogram()

# col_clusters_dat <- cutreeDynamic(col_hc,distM = as.matrix(col_mat),method="tree",minClusterSize = 12,
#                                   deepSplit = TRUE)
col_clusters_dat <- as.numeric(cutreeDynamic(col_hc,distM = as.matrix(col_mat),method="hybrid",
                                             minClusterSize = 9,
                                             deepSplit=1))
# col_clusters_dat <- cutree(col_hc,k=9)

# col_clusters_dat <- cutree(col_hc,k=5)
# col_clusters_dat8 <- cutree(col_hc,k=8)
col_clusters <- col_clusters_dat[order.dendrogram(col_dend)]
# mosaicplot(table(col_clusters,severity_data$model_cohort))


# plot(sapply(2:20,function(x){
#   out <- cutree(col_hc,k=x)
#   return(min(table(out)))
# }))

col_clusters_numbers <- unique(col_clusters) - (0 %in% col_clusters)
col_n_clusters <- length(col_clusters_numbers)
col_cols <- rainbow_hcl(col_n_clusters)
col_cols <- rep(brewer.pal(col_n_clusters,"Dark2"),times = table(col_cols))
# col_dend <- col_dend %>% branches_attr_by_clusters(col_clusters, values = col_cols)

# clustering protiens
row_mat <- severity_data[,pipeline_markers]
beta_k <- 1

# plot(as.numeric(cor(row_mat,method="pearson")),as.numeric(cor(row_mat,method="spearman")))
# abline(0,1)

# row_adj <- abs(cor(row_mat,method="pearson"))^beta_k
row_adj <- abs(cor(row_mat,method="spearman"))^beta_k

# row_adj <- cor(row_mat,method="pearson")^beta_k
# k <- as.vector(apply(row_adj,2,sum,na.rm=TRUE))
# par(mfrow=c(1,2))
# hist(k)
# scaleFreePlot(k,main="Check scale free topology")
# par(mfrow=c(1,1))

row_mat <- (1-row_adj) %>% as.dist()
row_hc <-  row_mat %>% hclust(method="ward.D2")
# row_hc <-  row_mat %>% hclust(method="complete")
row_dend <- row_hc %>% as.dendrogram()
# row_clusters_dat <- cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="tree",
#                                   minClusterSize = 10,deepSplit = TRUE)

# row_clusters_dat <- cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="tree",
#                                   minClusterSize = 10,deepSplit = TRUE)


# plot(sapply(2:20,function(x){
#   out <- cutree(row_hc,k=x)
#   return(min(table(out)))
# }))

# row_clusters_dat <- cutree(row_hc,k=6)
# row_clusters_dat <- cutree(row_hc,k=5)

# row_clusters_dat <- as.numeric(cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="tree",minClusterSize = 16))
row_clusters_dat <- as.numeric(cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="hybrid",minClusterSize = 15,
                                             deepSplit = TRUE))
# row_clusters_dat <- cutree(row_hc,k=4)


row_clusters <- row_clusters_dat[order.dendrogram(row_dend)]
row_clusters_numbers <- unique(row_clusters) - (0 %in% row_clusters)
row_n_clusters <- length(row_clusters_numbers)
row_cols <- rainbow_hcl(row_n_clusters)
row_cols <- rep(brewer.pal(row_n_clusters,"Dark2"),times = table(row_cols))
# row_dend <- row_dend %>% branches_attr_by_clusters(row_clusters, values = row_cols)

severity_heatmap <- t(severity_data[,pipeline_markers])
rownames(severity_heatmap) <- marker_to_gene(pipeline_markers)

# clean_clusters <-  factor(ifelse(row_clusters_dat==1,"CNS Regeneration",
#                                  ifelse(row_clusters_dat==2,"Related Innate Immunity","Adaptive Immunity")))
# module_colors <- c("CNS Regeneration" = "#1B9E77", "Related Innate Immunity" = "#D95F02","Adaptive Immunity"="#7570B3")
# clean_clusters <- row_clusters_dat
# module_colors <- c("1" = "#1B9E77", "2" = "#D95F02")
# draw(Heatmap(severity_heatmap,
#              cluster_columns = col_dend,cluster_rows = row_dend,
#              name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
#              show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
#              col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
#              column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm"),
#              row_names_gp = gpar(fontsize = 8),
#              top_annotation = HeatmapAnnotation("Cluster Solution 1" = as.factor(col_clusters_dat),
#                                                 "Cluster Solution 2" = as.factor(col_clusters_dat8))) +
#        rowAnnotation("Protein Clusters" = clean_clusters,
#                      col = list(`Protein Clusters` = module_colors)),
#      show_annotation_legend = FALSE)


severity <- tibble(soma = pipeline_markers)
severity_disc <- merge(severity,expression_dat,by="soma",all.x=TRUE)
severity_neuro <- severity_disc[,levs] %>% 
  as.matrix()
row.names(severity_neuro) <- severity_disc$gene

# col_clusters_dat <- as.numeric(cutreeDynamic(col_hc,distM = as.matrix(col_mat),method="tree",
#                                              minClusterSize = 9,
#                                              deepSplit=1))
# col_clusters_dat <- col_clusters_dat + 1



severity_first <- ifelse(pipeline_markers %in% unlist(strsplit(var_lists[[3]],"/")),"Yes","No")
severity_last <- ifelse(pipeline_markers %in% unlist(strsplit(var_lists[[4]],"/")),"Yes","No")
severity_first_mri <- ifelse(pipeline_markers %in% unlist(strsplit(var_lists[[7]],"/")),"Yes","No")

row_clusters_dat <- as.numeric(cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="hybrid",minClusterSize = 20,
                                             deepSplit = TRUE))
col_clusters_dat <- as.numeric(cutreeDynamic(col_hc,distM = as.matrix(col_mat),method="hybrid",
                                             minClusterSize = 12,
                                             deepSplit=1))
# col_clusters_dat <- cutree(col_hc,k=4)
table(col_clusters_dat)

# colss <- 1:length(unique(col_clusters_dat))
colss <- rep("white",length(unique(col_clusters_dat)))
names(colss) <- 1:length(unique(col_clusters_dat))
# colss <- list("Cluster" = colss)

marker_colss <- 1:length(unique(row_clusters_dat))
names(marker_colss) <- 1:length(unique(row_clusters_dat))
markers_cols <- list("Protein Clusters" = marker_colss)

severity_data <- severity_data %>% 
  mutate(diagnosis = str_replace_all(as.character(diagnosis),"-","")) %>% 
  mutate(diagnosis = ifelse(diagnosis =="HD","HV",diagnosis)) %>% 
  mutate(diagnosis = factor(diagnosis,levels = c("HV","RRMS","SPMS","PPMS")))


clus_cols <- list("Cluster" = colss,
                  Gender=c("Female"="pink","Male"="steelblue"),
                  Diagnosis = c("RRMS"="yellow","SPMS"="orange","PPMS"="brown"))



# col_dend2 <- col_dend %>% 
#   rotate(as.character(1:192))

row.names(severity_heatmap)[89] <- "PRKAG1"




tiff(my_filename("./results/severity_heatmap/severity_heatmap_V2.tiff"),units="in",res=300,width=14.5,height=10)
draw(severity_heatmap_obj <- Heatmap(severity_heatmap,
             cluster_columns = col_dend,
             cluster_rows = row_dend,
             # cluster_rows = FALSE,
             name="SomaScan expression (Z-score)",show_column_names = FALSE,
             show_row_names = TRUE,
             column_title = "Clustering of proteins that constitute three MS severity models",
             width=10, column_dend_reorder = FALSE,
             col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
             column_dend_height = unit(1, "cm"),row_dend_width = unit(1, "cm"),
             row_names_gp = gpar(fontsize = 6),row_names_side = "left",
             heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(5, "cm"), 
                                         title_position = "lefttop"),
             top_annotation = HeatmapAnnotation(data.frame("Cluster" = as.factor(col_clusters_dat),
                                                           "Gender" = severity_data$sex,
                                                           "Diagnosis" = severity_data$diagnosis),
                                                col = clus_cols,
                                                show_annotation_name = c(TRUE,TRUE,TRUE),
                                                annotation_name_side = c("left","left","left"))) +
       rowAnnotation("Protein Clusters" = as.factor(row_clusters_dat),
                     col=markers_cols) +
       rowAnnotation("MS-DSS: Baseline" = severity_first,
                     "MRI Severity: Baseline" = severity_first_mri,
                     "MS-DSS: Follow-up" = severity_last,
                     col = list("MS-DSS: Baseline"=c("Yes"="black","No"="white"),
                                "MRI Severity: Baseline"=c("Yes"="black","No"="white"),
                                "MS-DSS: Follow-up"=c("Yes"="black","No"="white")),
                     show_legend=FALSE,show_annotation_name=TRUE,annotation_name_side="top",
                     annotation_name_rot = 90),
     show_annotation_legend = TRUE,heatmap_legend_side="bottom",annotation_legend_side="right",
     padding = unit(c(2, 2, 6, 2),"mm"),gap=unit(c(2,15),"mm"))
dev.off()

col_dend2 <- column_dend(severity_heatmap_obj)[[1]] %>% 
  rotate(as.character(1:192))
plot(col_dend2)

names(row_dend(severity_heatmap_obj))
# patient_clusters <- tibble(patientcode = severity_data$patientcode,
#                            patient_clusters = c(9,7,5,3,6,4,8,1,2)[col_clusters_dat])
# write_csv(patient_clusters,
#           "./data/PatientClusters_20190524.csv")

severity_clusters <- tibble(marker = pipeline_markers,
                        cluster = row_clusters_dat)
severity_clusters <- severity_clusters[row_order(severity_heatmap_obj)[[1]],]

severity_clusters <- inner_join(clean_marker(severity_clusters$marker),severity_clusters,by="marker")
# write_csv(severity_clusters,my_filename("./results/string/pipeline_clusters/Pipeline_Markers_Clusters.csv"))

# cor_dat <- severity_data[,severity_clusters$marker]
# names(cor_dat) <- marker_to_gene(severity_clusters$marker)

# Bibi way
cor_col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061")))(200)

# cor_dat <- severity_data[,c("msdss","msdss_last","mri_severity",severity_clusters$marker)]
# cor_dat <- severity_data[,c("soma_msdss_adjust","soma_msdss_last_adjust","soma_mri_severity",severity_clusters$marker)]
cor_dat <- severity_data[,c("soma_msdss_adjust_v2","soma_msdss_last_adjust_v2","soma_mri_severity_v2",severity_clusters$marker)]
names(cor_dat) <- c("MS-DSS","MS-DSS (Followup)","MRI Severity",marker_to_gene(severity_clusters$marker))

cor_mat <- foreach(i = 1:3,.combine="cbind") %do% {
  sapply(marker_to_gene(severity_clusters$marker),function(x){as.numeric(cor.test(cor_dat[[i]],cor_dat[[x]],method="spearman")$est)})
}
cor_p <- foreach(i = 1:3,.combine="cbind") %do% {
  sapply(marker_to_gene(severity_clusters$marker),function(x){as.numeric(cor.test(cor_dat[[i]],cor_dat[[x]],method="spearman")$p.value)})
}

tiff("./results/severity_heatmap/please_work.tiff",width=1.5,height=10,units="in",res=300)
corrplot(cor_mat,is.corr=FALSE,method="square",tl.pos = "l",addgrid.col = NA,type="full",
         tl.col=c("black","red","green","blue")[c(as.numeric(as.factor(severity_clusters$cluster)))],
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(2,2,2,2),col=cor_col)
dev.off()


colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")
tiff("check.tiff",units="in",width=24,height=10,res=300)
draw(severity_heatmap_obj+Heatmap(cor(cor_dat[,-c(1:3)],method="spearman"),cluster_rows=FALSE,cluster_columns=FALSE,show_column_names = FALSE,
                                  row_names_side = "left",
             row_names_gp = gpar(fontsize = 5,col=c("black","red","green","blue")[c(as.numeric(as.factor(severity_clusters$cluster)))]),
             col=colorRamp2(c(seq(-1,1,length=200)),c(cor_col),space="sRGB"),width=6.5/13,name="left",
             heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(2, "cm"), title_position = "topleft"))+
       Heatmap(cor_mat,col = cor_col,cluster_columns = FALSE,cluster_rows = FALSE,row_dend_reorder = FALSE,row_names_side = "left",
               show_column_names = FALSE,show_row_names = FALSE,width=.75/13,name="right",
               heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(2, "cm"), title_position = "topleft",at=c(-0.4,0,.4))),
     heatmap_legend_side="bottom",
     show_annotation_legend = TRUE,annotation_legend_side="right",
     padding = unit(c(2, 2, 6, 2),"mm"))
dev.off()

tiff("check.tiff",units="in",width=8,height=8,res=300)
names(cor_dat)[99] <- "PRKAG1"
draw(Heatmap(cor(cor_dat[,-c(1:3)],method="spearman"),cluster_rows=FALSE,cluster_columns=FALSE,show_column_names = FALSE,
                                  row_names_side = "left",
                                  row_names_gp = gpar(fontsize = 6,col=c("black","red","green","blue")[c(as.numeric(as.factor(severity_clusters$cluster)))]),
                                  col=colorRamp2(c(seq(-1,1,length=200)),c(cor_col),space="sRGB"),width=6.5/8,name="left",
                                  heatmap_legend_param = list(
                                    legend_direction = "horizontal",legend_width = unit(5.9, "in"), title_position = "topleft"
                                    ))+
       Heatmap(cor_mat,col = cor_col,cluster_columns = FALSE,cluster_rows = FALSE,row_dend_reorder = FALSE,row_names_side = "left",
               show_column_names = FALSE,show_row_names = FALSE,width=1.5/8,name="right",
               heatmap_legend_param = list(
                 legend_direction = "horizontal",legend_width = unit(1.33, "in"), title_position = "topleft",at=c(-0.4,0,.4)
                 )),
     heatmap_legend_side="bottom",
     show_annotation_legend = TRUE,annotation_legend_side="right",
     padding = unit(c(2, 2, 6, 2),"mm"))
dev.off()


cor_dat <- severity_data[,pipeline_markers[row_order(severity_heatmap_obj)[[1]]]]
names(cor_dat) <- marker_to_gene(pipeline_markers[row_order(severity_heatmap_obj)[[1]]])

# names(cor_dat)[99] <- "PRKAG1"
names(cor_dat)[96] <- "PRKAG1"
rownames(cor_mat)[96] <- "PRKAG1"

library(ComplexHeatmap)
outcome_levs <- c("MS-DSS - Baseline","MRI Severity - Baseline", "MS-DSS - Follow-up")
tiff("./results/severity_heatmap/correlation_with_outcomes.tiff",width=1.5,height=10,units="in",res=300)
draw(Heatmap(cor_mat,col = cor_col,cluster_columns = FALSE,cluster_rows = FALSE,row_dend_reorder = FALSE,row_names_side = "left",
        show_column_names = FALSE,name="",
        row_names_gp = gpar(fontsize = 4,col=c("black","red","green","blue")[c(as.numeric(as.factor(severity_clusters$cluster)))]),
        # top_annotation = HeatmapAnnotation(
        #   text = anno_text(outcome_levs,rot=45,offset = unit(0, "npc"), just = "left",gp=gpar(fontsize=6)),
        #   annotation_height = max_text_width(outcome_levs)
        # ),
     heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(2, "cm"), title_position = "topleft",at=c(-0.4,0,.4))),
     heatmap_legend_side="bottom")
dev.off()

middle_64 <- cor_col[c(69:132)]

library(ComplexHeatmap)
outcome_levs <- c("MS-DSS - Baseline","MRI Severity - Baseline", "MS-DSS - Follow-up")
tiff("./results/severity_heatmap/correlation_with_outcomes_v2.tiff",width=1.5,height=10,units="in",res=300)
draw(Heatmap(cor_mat,col = cor_col[c(69:132)],cluster_columns = FALSE,cluster_rows = FALSE,row_dend_reorder = FALSE,row_names_side = "left",
             show_column_names = FALSE,name="Spearman",show_heatmap_legend = FALSE,
             row_names_gp = gpar(fontsize = 4,col=c("black","red","green","blue")[c(as.numeric(as.factor(severity_clusters$cluster)))]),
             # top_annotation = HeatmapAnnotation(
             #   text = anno_text(outcome_levs,rot=45,offset = unit(0, "npc"), just = "left",gp=gpar(fontsize=6)),
             #   annotation_height = max_text_width(outcome_levs)
             # ),
             heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(2, "cm"), title_position = "lefttop",at=c(-0.4,0,.4))),
     heatmap_legend_side="bottom")
dev.off()


# tiff(my_filename("./results/severity_heatmap/correlations_V2.tiff"),units="in",res=300,width=10,height=10)
# corrplot(cor(cor_dat,method="spearman",use="pairwise"),method="square",tl.pos = "l",addgrid.col = NA,type="full",
#          tl.col=c("black","red","green","blue","orange")[c(5,5,5,as.numeric(as.factor(severity_clusters$cluster)))],
#          tl.cex=0.45,title="Correlation among protein clusters",
#          mar = c(0, 0, 2, 0),col=cor_col)
# dev.off()

tiff(my_filename("./results/severity_heatmap/correlations_V3.tiff"),units="in",res=300,width=10,height=10)
corrplot(cor(cor_dat,method="spearman",use="pairwise"),method="square",tl.pos = "l",addgrid.col = NA,type="full",
         tl.col=c("black","red","green","blue")[c(as.numeric(as.factor(severity_clusters$cluster)))],
         tl.cex=0.45,title="Correlation among protein modules",
         mar = c(0, 0, 2, 0),col=cor_col,cl.pos = "b")
dev.off()


tiff(my_filename("./results/severity_heatmap/correlations_V3.tiff"),units="in",res=300,width=10,height=10)
corrplot(cor(cor_dat,method="spearman",use="pairwise"),method="square",tl.pos = "l",addgrid.col = NA,type="full",
         tl.col=c("black","red","green","blue")[c(as.numeric(as.factor(severity_clusters$cluster)))],
         tl.cex=0.45,title="Correlation among protein modules",
         mar = c(0, 0, 2, 0),col=cor_col,cl.pos = "n")
colorlegend()
dev.off()

plot(0,xlim = c(0,6), ylim = c(-0.5,1.2), type = "n")

colorlegend(cor_col,vertical=FALSE,labels=c(-.4,0,.4))

par(mfrow=c(2,2))
the_order <- "AOE"
corrplot(cor(severity_data[,c("mri_severity",var_lists[[7]])],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
corrplot(cor(severity_data[,c("mri_severity",unique(unlist(strsplit(var_lists[[7]],"/"))))],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
corrplot(cor(severity_data[,c("soma_mri_severity",var_lists[[7]])],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
corrplot(cor(severity_data[,c("soma_mri_severity",unique(unlist(strsplit(var_lists[[7]],"/"))))],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
par(mfrow=c(1,1))


par(mfrow=c(2,2))
the_order <- "original"
corrplot(cor(severity_data[,c("msdss",var_lists[[3]])],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
corrplot(cor(severity_data[,c("msdss",unique(unlist(strsplit(var_lists[[3]],"/"))))],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
corrplot(cor(severity_data[,c("soma_msdss_adjust",var_lists[[3]])],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
corrplot(cor(severity_data[,c("soma_msdss_adjust",unique(unlist(strsplit(var_lists[[3]],"/"))))],method="spearman",use="pairwise"),
         method="square",tl.pos = "lt",addgrid.col = NA,type="full",
         tl.cex=0.45,title="Correlation among protein clusters",
         mar = c(0, 0, 2, 0),col=cor_col,order=the_order)
par(mfrow=c(1,1))


clus1_pca <- princomp(severity_data[severity_clusters[severity_clusters$cluster == 1,][["marker"]]],cor=TRUE)
# soma_untreated_adj$module1 <- predict(clus1_pca,soma_untreated_adj)[,1]
module1 <- predict(clus1_pca,severity_data)[,1]

clus2_pca <- princomp(severity_data[severity_clusters[severity_clusters$cluster == 2,][["marker"]]],cor=TRUE)
# soma_untreated_adj$module2 <- predict(clus2_pca,soma_untreated_adj)[,1]
module2 <- predict(clus2_pca,severity_data)[,1]

clus3_pca <- princomp(severity_data[severity_clusters[severity_clusters$cluster == 3,][["marker"]]],cor=TRUE)
# soma_untreated_adj$module3 <- predict(clus3_pca,soma_untreated_adj)[,1]
module3 <- predict(clus3_pca,severity_data)[,1]

library(ggplot2)
fun_dat <- tibble(id = severity_data$patient,
                  cluster = col_clusters_dat,
                  module1=module1,
                  module2=module2,
                  module3=module3,
                  msdss = severity_data$msdss_last) %>% 
  gather(key,value,module1:module3) %>%
  # gather(key,value,module1:msdss) %>%
  group_by(key) %>% 
  mutate(value = (value - min(value))/(max(value)-min(value)))

fun_dat %>% 
  ggplot(aes(x=as.factor(cluster),y=value,color=key)) +
  geom_boxplot(notch = FALSE) +
  # geom_jitter(height = 0,width=0.05) +
  xlab("Patient Clusters") +
  ylab("Module Values (Scaled)") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")

Heatmap(t(soma_untreated_adj[soma_untreated_adj$diagnosis != "HD",pipeline_markers]),
        clustering_method_columns = "ward.D2",cluster_columns = TRUE,
        clustering_distance_columns = "euclidean",name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
        show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
        clustering_distance_rows = "pearson",clustering_method_rows = "average",
        col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
        column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm"))

# draw(Heatmap(t(soma_untreated_adj[,pipeline_markers]),
#         clustering_method_columns = "ward.D2",cluster_columns = TRUE,
#         clustering_distance_columns = "euclidean",name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
#         show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
#         clustering_distance_rows = "pearson",clustering_method_rows = "ward.D2",
#         col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
#         column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm"),
#         row_names_gp = gpar(fontsize = 8),
#         top_annotation = diag_at),annotation_legend_side="top")
# 
# draw(Heatmap(t(soma_untreated_adj[,pipeline_markers]),
#              clustering_method_columns = "ward.D2",cluster_columns = TRUE,
#              clustering_distance_columns = "euclidean",name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
#              show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
#              clustering_distance_rows = "pearson",clustering_method_rows = "ward.D2",
#              col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
#              column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm"),
#              row_names_gp = gpar(fontsize = 8),
#              top_annotation = diag_at),annotation_legend_side="right")

Heatmap(expression_dat[expression_dat$soma %in% pipeline_markers,-c(1:4)],
        clustering_method_columns = "ward.D2",cluster_columns = TRUE,
        clustering_distance_columns = "euclidean",name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
        show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
        clustering_distance_rows = "pearson",clustering_method_rows = "average",
        col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"))

##########################################################################
#################### Disability ##########################################
##########################################################################
disability <- read_excel("./data/Data for String and Zhang.xlsx",sheet = "disability")

# check <- filter(disability,soma %in% with(disability,soma[which(duplicated(soma))]))

disability <- disability %>% 
  mutate(source = ifelse(soma %in% with(disability,soma[which(duplicated(soma))]),"both",source)) %>% 
  remove_dups()
disability_disc <- merge(disability,expression_dat,by="soma")
disability_dir <- disability_disc %>% 
  filter(soma %in% disability_disc$soma) %>% 
  # mutate(direction = ifelse(direction=="negative","negative",paste(source,direction,sep="-"))) %>% 
  .$direction

disability_neuro <- disability_disc[,levs] %>% 
  as.matrix()

soma_untreated_adj <- soma_untreated_adj %>% 
  arrange(brainparfr1_volume)

disability_soma <- t(soma_untreated_adj[,disability_disc$soma])
row.names(disability_soma) <- disability_disc$gene
row.names(disability_neuro) <- disability_disc$gene

draw(Heatmap(disability_neuro,clustering_method_rows = "ward.D2",cluster_columns = TRUE,
             clustering_distance_rows = "euclidean",name="CNS Expression (Z-scores)",show_column_names = TRUE,
             show_row_names = TRUE,column_title = "CNS Expression",width = 1,
             col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),row_names_side = "left",
             split = disability_dir,row_names_gp = gpar(fontsize = 4)) +
       Heatmap(disability_soma,clustering_method_columns = "ward.D2",cluster_columns = TRUE,
               clustering_distance_columns = "euclidean",name="Somalogic Expression (Z-scores)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "MS Patient Expression",width=10,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
     column_title = "Proteins with Disability Associations",
     show_heatmap_legend=TRUE)

Heatmap(disability_soma,clustering_method_columns = "ward.D2",cluster_columns = TRUE,
        clustering_distance_columns = "euclidean",name="SomaScan Expression (Z-score)",show_column_names = TRUE,
        show_row_names = FALSE,column_title = "MS Patient Expression",width=10,
        clustering_distance_rows = "pearson",clustering_method_rows="ward.D2",
        col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),split=disability_dir)

#########################################################################
############# Sex Associations - Adjusted ###############################
#########################################################################
sex <- read_excel("./data/Data for String and Zhang.xlsx",sheet = "sex")
sex_disc <- merge(sex,expression_dat,by="soma")
sex_dir <- sex_disc %>% 
  filter(soma %in% sex_disc$soma) %>% 
  .$direction

sex_neuro <- sex_disc[,levs] %>% 
  as.matrix()

sex_soma <- t(soma_untreated_adj[,sex_disc$soma])
row.names(sex_soma) <- sex_disc$gene
row.names(sex_neuro) <- sex_disc$gene

sex_ha <- rowAnnotation(df = data.frame(Direction=sex_dir),
                        col = list(Direction = c("Elevated in males" = "blue", "Elevated in females" = "red")))

draw(Heatmap(sex_neuro,clustering_method_rows = "ward.D2",cluster_columns = TRUE,
             clustering_distance_rows = "euclidean",name="CNS Expression (Z-scores)",show_column_names = TRUE,
             show_row_names = TRUE,column_title = "CNS Expression",width = 1,
             col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),row_names_side = "left",
             split = sex_dir) +
       Heatmap(sex_soma[,soma_untreated_adj$sex=="Male"],clustering_method_columns = "ward.D2",cluster_columns = TRUE,
               clustering_distance_columns = "euclidean",name="Somalogic Expression (Z-scores)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "Male Patient Expression",width=4.5,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")) + 
       Heatmap(sex_soma[,soma_untreated_adj$sex=="Female"],clustering_method_columns = "ward.D2",cluster_columns = TRUE,
               clustering_distance_columns = "euclidean",name="Somalogic Expression (Z-scores)",show_column_names = TRUE,
               show_row_names = FALSE,column_title = "Female Patient Expression",width=5.5,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB")),
     column_title = "Proteins with Sex Associations in Healthy Donor Data",
     show_heatmap_legend=TRUE)


#################################################################
################### All somalogic - Unadjusted ##################
#################################################################
draw(Heatmap(expression_dat[,levs],clustering_method_rows = "ward.D2",cluster_columns = TRUE,
             clustering_distance_rows = "euclidean",name="Neuronal Expression of Zhang et al.",
             show_column_names = TRUE,show_row_names = TRUE,
             col=colorRamp2(c(-1.7,0,3.95),c("green","black","red"),space="sRGB"),width = 1) +
       Heatmap(t(soma_untreated_adj[,expression_dat$soma]),clustering_method_columns = "ward.D2",
               cluster_columns = TRUE,clustering_distance_columns = "euclidean",name="Somalogic Expression",
               show_column_names = TRUE,show_row_names = FALSE,column_title = "MS Patient Expression Data",
               width=5,
               col=colorRamp2(c(-7,0,8),c("green","black","red"),space="sRGB")),
     column_title = "All Somalogic Data")

# Clustering patients
col_mat <- soma_untreated[soma_untreated$diagnosis!="HD",expression_dat$soma] %>% dist
col_hc <-  col_mat %>% hclust(method="ward.D2")
# col_hc <-  col_mat %>% hclust(method="average")
col_dend <- col_hc %>% as.dendrogram()
col_clusters_dat <- cutreeDynamic(col_hc,distM = as.matrix(col_mat),method="hybrid",minClusterSize = 10)
# col_clusters_dat <- as.numeric(cutreeDynamic(col_hc,distM = as.matrix(col_mat),method="hybrid",minClusterSize = 10))
col_clusters <- col_clusters_dat[order.dendrogram(col_dend)]
col_clusters_numbers <- unique(col_clusters) - (0 %in% col_clusters)
col_n_clusters <- length(col_clusters_numbers)
col_cols <- rainbow_hcl(col_n_clusters)
col_cols <- rep(brewer.pal(col_n_clusters,"Dark2"),times = table(col_cols))
col_dend <- col_dend %>% branches_attr_by_clusters(col_clusters, values = col_cols)

# clustering protiens
row_mat <- soma_untreated[soma_untreated$diagnosis!="HD",expression_dat$soma]
beta_k <- 1
row_adj <- abs(cor(row_mat,method="pearson"))^beta_k
k <- as.vector(apply(row_adj,2,sum,na.rm=TRUE))
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check scale free topology")
par(mfrow=c(1,1))

row_mat <- (1-row_adj) %>% as.dist()
row_hc <-  row_mat %>% hclust(method="ward.D2")
# row_hc <-  row_mat %>% hclust(method="average")
row_dend <- row_hc %>% as.dendrogram()
row_clusters_dat <- cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="tree",
                                  minClusterSize = 20,deepSplit = TRUE)

row_clusters <- row_clusters_dat[order.dendrogram(row_dend)]
row_clusters_numbers <- unique(row_clusters) - (0 %in% row_clusters)
row_n_clusters <- length(row_clusters_numbers)
row_cols <- rainbow_hcl(row_n_clusters)
# row_cols <- rep(brewer.pal(row_n_clusters,"YlGnBu"),times = table(row_cols))
row_dend <- row_dend %>% branches_attr_by_clusters(row_clusters, values = row_cols)

library(WGCNA)
row_clusters_dat <- cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="tree",
                                  minClusterSize = 20,deepSplit = FALSE)
row_clusters_dat2 <- cutreeDynamic(row_hc,distM = as.matrix(row_mat),method="tree",
                                   minClusterSize = 20,deepSplit = TRUE)
row_cols <- rep(brewer.pal(11,"Spectral"),times = table(row_clusters_dat))
plotDendroAndColors(row_hc, colors = data.frame(module = labels2colors(row_clusters_dat),
                                                deep=labels2colors(row_clusters_dat2)),
                    dendroLabels = FALSE, 
                    main = "Gene dendrogram and module colors")
table(row_clusters_dat)
table(row_clusters_dat2)

Heatmap(t(soma_untreated[soma_untreated$diagnosis!="HD",expression_dat$soma]),
               cluster_columns = col_dend,cluster_rows=row_dend,name="Somalogic Expression",
               show_column_names = FALSE,show_row_names = FALSE,column_title = "MS Patient Expression Data",
               width=5,
               col=colorRamp2(c(-7,0,8),c("green","black","red"),space="sRGB"),
               column_dend_height = unit(3, "cm"),row_dend_width = unit(3, "cm"))

soma_clusters <- tibble(marker = expression_dat$soma,
                        cluster = row_clusters_dat) %>% 
  arrange(cluster)


library(corrplot)
corrplot(cor(soma_untreated[,soma_clusters$marker]),method="square",tl.pos = "n",addgrid.col = NA)


clinical_outcomes <- c("armss_last","msdss_last","armss","msdss",
                       "combiwise_score","neurex_total","comris_ctd")
volumes <- c("mri_severity","volume_calculated_brainparenchymalfr","volume_calculated_csf_fr",
             "volume_calculated_gmfr","volume_calculated_wmfr","volume_calculated_lesionfr",
             "volume_calculated_nawmfr",
             "volume_calculated_supratentorialgm","volume_calculated_supratentorialwm",
             "volume_calculated_totalbrain","volume_calculated_totalcsf")
outcomes_to_check <- c("diagnosis","gender","age",clinical_outcomes[-c(1)],volumes,"csf_protein_total","csf_iggindex",
                       "csf_albuminquotient",qual)

outcomes <- c("age",clinical_outcomes,volumes)
# outcomes <- c("age",clinical_outcomes)

rank_inverse <- function(x,k=3/8){
  n <- sum(!is.na(x))
  out <- (rank(x)-k)/(n-2*k+1)
  return(qnorm(out))
}
my_scale <- function(x){
  x <- (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
  return(x)
}
pat_clust_outcome <- soma_untreated_adj %>% 
  filter(diagnosis !="HD") %>% 
  filter(sampleid != "NDUS1187") %>% 
  mutate(patient_clusters = col_clusters_dat) %>% 
  # mutate_at(outcomes,my_scale)
  # mutate_at(c(clinical_outcomes,volumes),z_score,log=TRUE)
  # mutate_at(outcomes,z_score,log=FALSE)
  mutate_at(outcomes,rank_inverse) %>% 
  mutate_at(outcomes,my_scale)

library(cluster)
col_dend <- as.dendrogram(hclust(daisy(pat_clust_outcome[,c(outcomes,qual)],metric="gower"),method="ward.D2"))

pat_clust_outcome <- pat_clust_outcome %>% 
  mutate_at(qual,as.numeric) %>%
  mutate_at(qual,my_scale)

row_dend <- as.dendrogram(hclust(as.dist(1-abs(cor(pat_clust_outcome[,outcomes],method="pearson",use = "pairwise"))),method="ward.D2"))
# row_dend <- as.dendrogram(hclust(as.dist(1-cor(pat_clust_outcome[,outcomes],method="pearson",use = "pairwise")),method="ward.D2"))
library(StatMatch)
clean_dist <- function(x,y){
  if(mean(class(x)==class(y))!=1){stop("Variable classes don't match!")}
  if(class(x)=="numeric"){1-abs(cor(x,y,use="pairwise"))}
  else{
    if("Mild" %in% unique(x)){x <- ordered(x,levels=c("None","Mild","Moderate","Severe"))}
    if("Mild" %in% unique(y)){y <- ordered(y,levels=c("None","Mild","Moderate","Severe"))}
    if("<50%" %in% unique(x)){x <- ordered(x,levels=c("None","<50%",">50%"))}
    if("<50%" %in% unique(y)){y <- ordered(y,levels=c("None","<50%",">50%"))}
    # print(x)
    # print(y)
    x <- as.numeric(x)
    y <- as.numeric(y)
    return(gower.dist(x,y,KR.corr = FALSE))}
}
full_outcomes <- c(outcomes,qual)
to_split <- ifelse(full_outcomes %in% outcomes,"Quantitative","Semi-Quantiative")
qual_red <- c("lesion_load_brainstem","lesion_load_cerebellum","lesion_load_medulla")
ht_list <- sapply(sort(unique(col_clusters_dat)),function(i){
  add_r <- FALSE
  if(i == 9){add_r <- TRUE}
  Heatmap(t(pat_clust_outcome[col_clusters_dat==i,c(outcomes,qual_red)]),
          cluster_columns = TRUE,
          row_names_side = "right",
          show_row_names = add_r,
          show_column_names = FALSE,
          cluster_rows=TRUE,
          na_col = "orange",
          name="Standardized Outcomes",
          clustering_distance_rows = "pearson",clustering_method_rows = "ward.D2",
          clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D2",
          col=colorRamp2(c(0,0.5,1),c("blue","black","red"),space="sRGB"),
          column_dend_height = unit(2, "cm"),heatmap_legend_side = "bottom",
          column_dend_reorder = TRUE,row_dend_reorder = TRUE,
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(direction = "horizontal"),
          top_annotation = HeatmapAnnotation(data.frame(Sex = pat_clust_outcome$sex[col_clusters_dat==i],
                                                        Diagnosis = pat_clust_outcome$diagnosis[col_clusters_dat==i]),
                                             col = list(Sex=c("Female"="green","Male"="purple"),
                                                        Diagnosis = c("RR-MS"="red","SP-MS"="green","PP-MS"="yellow")),
                                             show_annotation_name=add_r,annotation_name_side="right"))
},simplify = FALSE)
draw(Reduce(`+`,ht_list),heatmap_legend_side = "bottom",
     column_title = "Examining outcomes across patient clusters (Limited qualitative MRI)")


Heatmap(t(pat_clust_outcome[,c(outcomes,qual)]),
        cluster_columns = TRUE,
        row_names_side = "right",
        cluster_rows=TRUE,
        split=to_split,
        clustering_distance_rows = "pearson",clustering_method_rows = "ward.D2",
        clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D2",
        col=colorRamp2(c(0,0.5,1),c("green","black","red"),space="sRGB"),
        column_dend_height = unit(2, "cm"),
        column_dend_reorder = TRUE,row_dend_reorder = TRUE,
        show_heatmap_legend = TRUE)



# code only for demonstration
ht_list = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
for(s in sth) {
  ht_list = ht_list + Heatmap(...)
}


Heatmap(t(pat_clust_outcome[,c(outcomes,qual)]),
        # cluster_columns = column_dend(severity_heatmap_obj)[[1]],
        cluster_columns = TRUE,
        row_names_side = "right",
        cluster_rows=TRUE,
        split=to_split,
        clustering_distance_rows = "pearson",clustering_method_rows = "ward.D2",
        clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D2",
        # col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
        col=colorRamp2(c(0,0.5,1),c("green","black","red"),space="sRGB"),
        column_dend_height = unit(2, "cm"),
        column_dend_reorder = TRUE,row_dend_reorder = TRUE,
        top_annotation = HeatmapAnnotation("Cluster Solution" = pat_clust_outcome$patient_clusters,
                                           col = clus_cols,show_legend = FALSE),
        show_heatmap_legend = TRUE)

Heatmap(t(pat_clust_outcome[,c(qual)]),
        # cluster_columns = column_dend(severity_heatmap_obj)[[1]],
        cluster_columns = col_dend,
        row_names_side = "right",
        cluster_rows=TRUE,
        # split=to_split,
        clustering_distance_rows = clean_dist,clustering_method_rows = "ward.D2",
        clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D2",
        col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
        column_dend_height = unit(2, "cm"),
        column_dend_reorder = TRUE,row_dend_reorder = TRUE,
        top_annotation = HeatmapAnnotation("Cluster Solution" = pat_clust_outcome$patient_clusters,
                                           col = clus_cols,show_legend = FALSE),
        show_heatmap_legend = TRUE)


row_dend <- as.dendrogram(hclust(as.dist(1-abs(cor(pat_clust_outcome[,outcomes],method="spearman"))),method="ward.D2"))

# Heatmap(t(pat_clust_outcome[,outcomes]),
#         cluster_columns = col_dend,cluster_rows=row_dend,
#         column_dend_height = unit(4, "cm"),
#         column_dend_reorder = FALSE,row_dend_reorder = TRUE,
#         clustering_distance_rows = "spearman",clustering_method_rows = "ward.D2",
#         col=colorRamp2(c(0,0.5,1),c("green","black","red"),space="sRGB"),
#         # top_annotation = HeatmapAnnotation("Cluster Solution" = as.factor(col_clusters_dat),
#         top_annotation = HeatmapAnnotation("Cluster Solution" = pat_clust_outcome$patient_clusters,
#                                            col = clus_cols,show_legend = FALSE))
Heatmap(t(pat_clust_outcome[,outcomes]),
        cluster_columns = column_dend(severity_heatmap_obj)[[1]],cluster_rows=row_dend,
        column_dend_height = unit(2, "cm"),
        column_dend_reorder = FALSE,row_dend_reorder = FALSE,
        clustering_distance_rows = "spearman",clustering_method_rows = "ward.D2",
        col=colorRamp2(c(0,0.5,1),c("green","black","red"),space="sRGB"),
        # top_annotation = HeatmapAnnotation("Cluster Solution" = as.factor(col_clusters_dat),
        top_annotation = HeatmapAnnotation("Cluster Solution" = pat_clust_outcome$patient_clusters,
                                           col = clus_cols,show_legend = FALSE),
        show_heatmap_legend = FALSE)


Heatmap(pat_clust_outcome[,outcomes],
        # cluster_columns = col_dend,
        clustering_distance_columns = "pearson",clustering_method_columns = "ward.D2",
        clustering_distance_rows = "euclidean",clustering_method_rows = "ward.D2",
        col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
        split = as.factor(col_clusters_dat),width=5)


qual

Heatmap(t(pat_clust_outcome[,qual]),
        cluster_columns = column_dend(severity_heatmap_obj)[[1]],
        row_names_side = "right",cluster_rows=TRUE,
        clustering_distance_rows = "pearson",clustering_method_rows = "ward.D2",
        clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D2",
        col=colorRamp2(c(0,0.5,1),c("green","black","red"),space="sRGB"),
        # col = brewer.pal(6,"RdGy"),
        column_dend_height = unit(2, "cm"),
        # column_dend_reorder = TRUE,row_dend_reorder = TRUE,
        column_dend_reorder = TRUE,row_dend_reorder = TRUE,
        # top_annotation = HeatmapAnnotation("Cluster Solution" = as.factor(col_clusters_dat),
        top_annotation = HeatmapAnnotation("Cluster Solution" = pat_clust_outcome$patient_clusters,
                                           col = clus_cols,show_legend = FALSE),
        show_heatmap_legend = TRUE)

patient_clusters <- tibble(patientcode = severity_data$patientcode,
                           patient_clusters = c(3,4,9,7,1,8,2,5,6)[col_clusters_dat])
# write_csv(patient_clusters,
#           "./data/PatientClusters_20190521.csv")

# soma_untreated <- soma_untreated %>% 
#   merge(patient_clusters,by="patientcode",sort = FALSE)
soma_untreated$patient_clusters <- as.factor(col_clusters_dat)


rank_inverse <- function(x,k=3/8){
  n <- sum(!is.na(x))
  out <- (rank(x)-k)/(n-2*k+1)
  return(qnorm(out))
}

hist(rgamma(1000,1,1))
hist(rank_inverse(rgamma(1000,1,1)))

# draw(Heatmap(severity_heatmap,
#              cluster_columns = col_dend,
#              cluster_rows = row_dend,
#              # cluster_rows = FALSE,
#              name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
#              show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
#              col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
#              column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm"),
#              row_names_gp = gpar(fontsize = 6),row_names_side = "left",
#              top_annotation = HeatmapAnnotation("Cluster Solution" = as.factor(col_clusters_dat),
#                                                 col = clus_cols)) +
#        rowAnnotation("Protein Clusters" = as.factor(clean_clusters),
#                      col=markers_cols),
#      show_annotation_legend = FALSE)


draw(Heatmap(severity_heatmap,
             cluster_columns = col_dend,
             cluster_rows = row_dend,
             # cluster_rows = FALSE,
             name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
             show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
             col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
             column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm"),
             row_names_gp = gpar(fontsize = 4),row_names_side = "left",
             top_annotation = HeatmapAnnotation("Cluster Solution" = as.factor(col_clusters_dat))) + 
       Heatmap(severity_neuro,clustering_method_rows = "ward.D2",cluster_columns = TRUE,
               clustering_distance_rows = "euclidean",name="CNS Expression (Z-scores)",show_column_names = TRUE,
               show_row_names = TRUE,column_title = "CNS Expression",width = 1,
               col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),row_names_side = "right",
               row_names_gp = gpar(fontsize = 4),
               column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm")) +
       rowAnnotation("Protein Clusters" = as.factor(clean_clusters)),
     show_annotation_legend = TRUE)

# draw(Heatmap(severity_neuro,clustering_method_rows = "ward.D2",cluster_columns = TRUE,
#              clustering_distance_rows = "euclidean",name="CNS Expression (Z-scores)",show_column_names = TRUE,
#              show_row_names = TRUE,column_title = "CNS Expression",width = 1,
#              col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),row_names_side = "left",
#              row_names_gp = gpar(fontsize = 4),
#              column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm")) + 
#        Heatmap(severity_heatmap,
#              cluster_columns = col_dend,cluster_rows = row_dend,
#              name="Somalogic Expression (Z-scores)",show_column_names = FALSE,
#              show_row_names = TRUE,column_title = "Proteins with Severity Associations",width=10,
#              col=colorRamp2(c(-3,0,3),c("green","black","red"),space="sRGB"),
#              column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm"),
#              row_names_gp = gpar(fontsize = 8),
#              top_annotation = HeatmapAnnotation("Cluster Solution 1" = as.factor(col_clusters_dat),
#                                                 "Cluster Solution 2" = as.factor(col_clusters_dat8))) + 
#        rowAnnotation("Protein Clusters" = clean_clusters,
#                      col = list(`Protein Clusters` = module_colors)),
#      show_annotation_legend = FALSE)

patient_clusters <- tibble(patientcode = severity_data$patientcode,
                           patient_clusters = c(3,4,9,7,1,8,2,5,6)[col_clusters_dat])
# write_csv(patient_clusters,
#           "./data/PatientClusters_20190521.csv")

soma_untreated <- soma_untreated %>% 
  left_join(patient_clusters,by="patientcode")

soma_untreated$patient_clusters <- as.factor(soma_untreated$patient_clusters)

# write_csv(tibble(patientcode = severity_data$patientcode,patient_clusters = c(4,1,7,5,6,2,3)[col_clusters_dat]),
#           "./data/PatientClusters_20190516.csv")


# write_csv(tibble(patient = severity_data$patient,patient_clusters = col_clusters_dat),
#           "./data/PatientClusters_FourthStab_20190425.csv")

