source("./scripts/0-functions.R")
source("./scripts/0-package-funs.R")

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
translate_file <- read_excel("./data/translation1.3k.xlsx")
tomatch <- translate_file %>% 
  rename(gene = EntrezGeneSymbol,protein=Target,protein_full=TargetFullName,marker=SomaId,id=SeqId) %>% 
  mutate(id = paste(gene,".",id,sep="")) %>% 
  mutate(id = str_replace_all(id,"-",".")) %>% 
  mutate(id = str_replace_all(id,"_",".")) %>% 
  mutate(id = str_replace_all(id,",",".")) %>%
  mutate(id = str_replace_all(id,"@",".")) %>%
  mutate(id = str_replace_all(id," ",".")) %>%
  mutate(id = str_replace_all(id,"(\\.)(\\.)",".")) %>%
  select(-UniProt,-EntrezGeneID)

soma_untreated <- read_csv("./data/CLEANED_soma_untreated_20190508.csv")
soma_untreated <- soma_untreated %>% 
  mutate(gender = factor(gender,levels = c("M","F"))) %>%
  mutate(mri_date = mdy(mri_date),dob=mdy(dob)) %>% 
  mutate(age_mri = as.numeric(mri_date-dob)/365)

new_markers <- soma_untreated %>% 
  select(starts_with("SL")) %>% 
  names()

age_to_adjust <- readLines("./data/age_to_adjust_20190626.txt")
sex_to_adjust <- readLines("./data/sex_to_adjust_20190626.txt")

new <- read_excel("./data/HD Information from Peter - 20181109.xlsx",sheet="1.3k HD")
new <- new[which(!duplicated(with(new,interaction(PatientCode,LPDate)))),]
names(new)[1:10] <- tolower(names(new))[1:10]

new <- new %>% 
  mutate(gender = factor(gender,levels=c("M","F")))

soma_new_adj <- foreach(i=1:length(new_markers),.final=as.data.table,.packages=c("dplyr")) %do% {
  mod <- adjust_function(new_markers[i],data=new)
  (log(soma_untreated[[new_markers[i]]]) - predict(mod,newdata=soma_untreated))
}
names(soma_new_adj) <- new_markers
soma_new_adj[,sampleid := soma_untreated$sampleid]
setcolorder(soma_new_adj,c("sampleid"))
soma_new_adj <- as_tibble(soma_new_adj)

soma_untreated_adj <- inner_join(select(soma_untreated,-starts_with("SL")),soma_new_adj,by="sampleid")

mri_mod <- lm(volume_calculated_brainparenchymalfr~age_mri,data=soma_untreated_adj)
soma_untreated_adj$mri_severity <- -1 * (soma_untreated_adj$volume_calculated_brainparenchymalfr - predict(mri_mod,newdata=soma_untreated_adj))

# diffs  <- c("t1_black_hole_fraction","t2_lesion_load",
#             "atrophy_type","cel_lesion_load","patient","mri_date","cel_exact_number",
#             "atrophy_score","atrophy_score_divide_by_cerebrum",
#             "lesion_score","lesion_score_plus_cerebrum","gm_plus_juxtacortical")
# clean_levels <- function(x,factorize = TRUE){
#   out <- str_replace_all(unlist(lapply(strsplit(x,"[(]"),first))," ","")
#   if(factorize==TRUE){out <- factor(out,levels=c("None","Mild","Moderate","Severe"))}
#   return(out)
# }
# 
# qual <- c("t1_black_hole_fraction","t2_lesion_load","cel_lesion_load","juxtacortical_lesion_load",
#           "lesion_atrophy_brainstem","lesion_load_brainstem","gm_lesion_load",
#           "lesion_load_cerebellum","lesion_atrophy_cerebellum",
#           "lesion_atrophy_medulla",
#           "lesion_load_medulla",
#           "lesion_cerebrum")
# 
# soma_untreated_adj <- soma_untreated_adj %>% 
#   mutate_at(names(soma_untreated_adj)[names(soma_untreated_adj) %in% qual& names(soma_untreated_adj) %notin% diffs],clean_levels)
# soma_untreated_adj$t1_black_hole_fraction <- factor(soma_untreated_adj$t1_black_hole_fraction,levels=c("None","<50%",">50%")) 
# soma_untreated_adj$t2_lesion_load <- factor(clean_levels(soma_untreated_adj$t2_lesion_load,factorize = FALSE),
#                                         levels = c("None","Mild","Moderate","Severe"))
# # soma_untreated$atrophy_type <- factor(soma_untreated$atrophy_type)
# soma_untreated_adj$cel_lesion_load <- factor(clean_levels(soma_untreated_adj$cel_lesion_load,factorize = FALSE),
#                                          levels=c("None","Mild","Moderate","Severe"))

outcomes <- c("combiwise_score","msss","armss","msdss","comris_ctd")
outcomes <- c(outcomes,paste(outcomes,"last",sep="_"))
outcomes <- c(outcomes,c("combiwise_score_slope_unadj","combiwise_score_slope_adj","neurex_total_slope_unadj","neurex_total_slope_adj"))
outcomes <- outcomes[-c(10)]

volumes <- str_subset(names(soma_untreated_adj),"volume")
volumes <- names(which(apply(soma_untreated_adj[,volumes],2,function(x){sum(is.na(x))})==23))
volumes <- volumes[!str_detect(volumes,"source")]

outcomes <- c(outcomes,volumes,"mri_severity")

# loads <- str_subset(names(soma_untreated_adj),"load")
# loads <- loads[!str_detect(loads,"comment")]
# loads <- loads[!str_detect(loads,"score")]
# atrophy <- str_subset(names(soma_untreated_adj),"atrophy")
# atrophy <- atrophy[!str_detect(atrophy,"comment")]
# atrophy <- atrophy[!str_detect(atrophy,"score")]
# atrophy <- atrophy[!str_detect(atrophy,"type")]
# atrophy <- atrophy[!str_detect(atrophy,"severity")]
# qual <- c(loads,atrophy)

##############################################################################
######## Individual Corrleations between Somamers and outcomes ###############
##############################################################################

# 1.3K Platform, Age/sex adjusted somamers
new_adj_est <- sapply(outcomes,function(ychar){
  sapply(new_markers,function(xchar){
  cor_est(soma_untreated_adj[[xchar]],soma_untreated_adj[[ychar]],method="spearman")
    })
  })
new_adj_pval <- sapply(outcomes,function(ychar){sapply(new_markers,function(xchar){
  cor_pval(soma_untreated_adj[[xchar]],soma_untreated_adj[[ychar]],method="spearman")})})

new_adj_pval <- as_tibble(new_adj_pval) %>% 
  mutate_all(function(x){p.adjust(x,method="fdr")}) %>%
  mutate(marker = new_markers,type="p-value") %>% 
  select(marker,type,everything())

new_adj_est <- as_tibble(new_adj_est) %>% 
  mutate(marker = new_markers,type="spearman") %>% 
  select(marker,type,everything())

new_adj <- bind_rows(new_adj_est,new_adj_pval)
new_adj <- new_adj %>% 
  gather(key,value,outcomes)         

new_adj <- new_adj %>% 
  mutate(key = factor(key,levels=unique(new_adj$key))) %>% 
  mutate(type=factor(type,levels=c("spearman","p-value"))) %>% 
  arrange(marker,key,type) %>% 
  mutate(key = paste(as.character(key),type)) %>% 
  select(-type) 

new_adj <- new_adj %>% 
  mutate(key = factor(key,levels=unique(new_adj$key))) %>% 
  spread(key,value)

new_adj <- right_join(tomatch,new_adj,by="marker")
write_csv(new_adj,my_filename("./results/soma_adjusted_outcome_correlations.csv"))

# jpeg(my_filename("./results/individual_correlations/MRI_severity_calculation.JPEG"),
#      res=300,units="in",width=6,height=5)
# plot(volume_calculated_brainparenchymalfr ~ age_mri, data = soma_untreated_adj,
#      main="MRI Severity",xlab="Age",ylab="Brain Parenchymal Fraction")
# abline(mri_mod)
# dev.off()

# par(mfrow=c(1,2))
# with(new_adj,plot(`combiwise_score_slope_adj spearman`,`neurex_total_slope_adj spearman`,asp=1,xlim=c(-.3,.3),ylim=c(-.3,.3)))
# with(new_adj,plot(`combiwise_score_slope_adj p-value`,`neurex_total_slope_adj p-value`,asp=1))
# par(mfrow=c(1,1))
# qual_pval_prop <- sapply(qual,function(ychar){sapply(new_markers,function(xchar){
#   kruskal.test(soma_untreated_adj[[xchar]]~soma_untreated_adj[[ychar]])$p.value})})
# 
# qual_pval_inaprop <- sapply(qual,function(ychar){sapply(new_markers,function(xchar){
#   cor_pval(soma_untreated_adj[[xchar]],as.numeric(soma_untreated_adj[[ychar]]),method="spearman")})})
# 
# qual_pval_prop <- as_tibble(qual_pval_prop) %>% 
#   mutate_all(function(x){p.adjust(x,method="fdr")}) %>%
#   mutate(marker = new_markers,type="Kruskal") %>% 
#   select(marker,type,everything())
# qual_pval_inaprop <- as_tibble(qual_pval_inaprop) %>% 
#   mutate_all(function(x){p.adjust(x,method="fdr")}) %>%
#   mutate(marker = new_markers,type="Spearman") %>% 
#   select(marker,type,everything())
# 
# qual_pval <- bind_rows(qual_pval_prop,qual_pval_inaprop)
# qual_pval <- qual_pval %>% 
#   gather(key,value,qual) %>% 
#   mutate(key = factor(key,levels = qual),
#          type = factor(type,levels=c("Kruskal","Spearman"))) %>% 
#   arrange(marker,key,type) %>% 
#   spread(type,value) %>% 
#   mutate(Kruskal_deteced = ifelse(Kruskal <= 0.05,"Yes","No"),
#          spearman_detected = ifelse(Spearman <= 0.05,"Yes","No"))
# 
# write_csv(qual_pval,"letsexplore.csv")
