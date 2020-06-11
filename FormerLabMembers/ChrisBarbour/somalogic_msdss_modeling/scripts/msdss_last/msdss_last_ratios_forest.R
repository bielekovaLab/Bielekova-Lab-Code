source("./scripts/0-functions.R")
source("./scripts/0-package-funs.R")

data <- read_csv("./data/CLEANED_soma_untreated_20190508.csv")
soma_dt <- readRDS("./data/soma_adjusted_ratios_20190626.rds")

data <- data %>% 
  filter(model_cohort=="training") %>% 
  select(sampleid,msdss_last) %>% 
  as.data.table()

soma_dt <- merge(data,soma_dt,by="sampleid",sort=TRUE)
soma_dt <- soma_dt[,-1]

all_rats <- names(soma_dt)

# INPUTS GO HERE
remove_outliers <- "FALSE"
var_path <- "./scripts/msdss_last/iter0_ratios.txt"
parseCommandArgs()
if(var_path != ""){
	var_list <- readLines(var_path)
	soma_dt <- subset(soma_dt,select=which(all_rats %in% c("msdss_last",var_list)))
 }

# Selecting tuning parameter
soma_dt <- as.data.frame(soma_dt)
nvars <- floor(3*sqrt(dim(soma_dt)[2]))
if(nvars >= dim(soma_dt)[2]-1){nvars <- floor(sqrt(dim(soma_dt)[2]))}
hold <- DUMX

ncpus <- detectBatchCPUs()
mod <- ranger(data=soma_dt, num.trees=40000,mtry=nvars,importance='impurity',num.threads=ncpus,
	      seed=DUMZ,dependent.variable.name="msdss_last")

saveRDS(mod,"DUMY1")
