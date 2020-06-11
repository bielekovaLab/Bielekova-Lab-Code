source("./scripts/0-functions.R")
source("./scripts/0-package-funs.R")

age_to_adjust <- readLines("./data/age_to_adjust_20190626.txt")
sex_to_adjust <- readLines("./data/sex_to_adjust_20190626.txt")

soma_untreated <- read_csv("./data/CLEANED_soma_untreated_20190508.csv")
new_markers <- soma_untreated %>% 
  select(starts_with("SL")) %>% 
  names()

soma_new <- soma_untreated %>% 
  select(sampleid,age,gender,which(names(soma_untreated) %in% new_markers)) %>% 
  mutate(gender = ifelse(gender=="Male","M","F")) 

new <- read_excel("./data/HD Information from Peter - 20181109.xlsx",sheet="1.3k HD")
new <- new[which(!duplicated(with(new,interaction(PatientCode,LPDate)))),]
names(new)[1:10] <- tolower(names(new))[1:10]

new <- new %>% 
  mutate(gender = factor(gender,levels=c("M","F")))

ratios_new <- character(length=choose(length(new_markers),2))
iter <- 1
for(i in 1:(length(new_markers)-1)){
  for(j in (i+1):length(new_markers)){
    ratios_new[iter] <- paste(new_markers[i],"/",new_markers[j],sep="")
    iter <- iter + 1
  }
}

soma_new_adj <- foreach(i=1:length(new_markers),.final=as.data.table,.packages=c("dplyr")) %do% {
  mod <- adjust_function(new_markers[i],data=new)
  (log(soma_new[[new_markers[i]]]) - predict(mod,newdata=soma_new))
}
names(soma_new_adj) <- new_markers
soma_new_adj[,sampleid := soma_new$sampleid]
setcolorder(soma_new_adj,c("sampleid"))

soma_dt_new <- foreach(i=1:length(ratios_new),.final=as.data.table,.packages=c('dplyr')) %do% {
  create_rat(ratios_new[i],dat=soma_new_adj,adj=TRUE)
}
names(soma_dt_new) <- ratios_new

soma_dt_new[,sampleid := soma_new_adj$sampleid]
setcolorder(soma_dt_new,c("sampleid"))

soma_dt_new <- merge(soma_dt_new,soma_new_adj,by="sampleid")
saveRDS(soma_dt_new,"./data/soma_adjusted_ratios_20190626.rds")

writeLines(c(ratios_new,new_markers),"./data/ratios_and_markers.txt")