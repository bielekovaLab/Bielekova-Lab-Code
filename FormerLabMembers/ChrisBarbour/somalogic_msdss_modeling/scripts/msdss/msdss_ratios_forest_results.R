source("./scripts/0-functions.R")
source("./scripts/0-package-funs.R")

data <- read_csv("./data/CLEANED_soma_untreated_20190508.csv")

data <- data %>% 
  filter(model_cohort=="training") %>% 
  select(sampleid,msdss) %>% 
  arrange(sampleid)

# INPUTS POTENTIALLY GOIONG TO COMMAND LINE GO HERE
remove_outliers <- "FALSE"
path1 <- "./scripts/msdss/msdss_ratios_forest_iter0_"
pull_iter <- "./scripts/msdss/iter0"
push_iter <- "./scripts/msdss/iter1"
parseCommandArgs()
y <- data$msdss

grab_oob <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(mod$prediction.error)
}

grab_cor <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(cor(y,mod$predictions))
}
ccc_function <- function(x1,x2){
	n <- length(x1)
	xb <- mean(x1)
	yb <- mean(x2)
	sx <- mean((x1-xb)^2)
	sy <- mean((x2-yb)^2)
	sxy <- mean((x1-xb)*(x2-yb))
	out <- (2*sxy)/(sx + sy + (xb - yb)^2)
	return(out)
}
grab_ccc <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(ccc_function(y,mod$predictions))
}

grab_imp <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(mod$variable.importance)
}

grab_imp_zero <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(mean(mod$variable.importance==0))
}
top_n <- function(x,n=500){
	 len<-length(x)
	 if(len < n){stop("n must be less than or equal to length of vector.")}
	 return(ifelse(x>quantile(x,1-n/len),"yes","no"))
}

pick_top_n <- function(x,n=500){
	 return(which(top_n(x,n)=="yes"))
}

path2 <- ".rds"
mod_list <- 1:10
prop_top <- 0.9

split_oob <- sapply(mod_list,grab_oob)
write.table(split_oob,paste(pull_iter,"_oob_results.txt",sep=""),row.names=FALSE)

split_cor <- sapply(mod_list,grab_cor)
write.table(split_cor,paste(pull_iter,"_oob_cor_results.txt",sep=""),row.names=FALSE)

split_ccc <- sapply(mod_list,grab_ccc)
write.table(split_ccc,paste(pull_iter,"_oob_ccc_results.txt",sep=""),row.names=FALSE)

split_imp <- sapply(mod_list,grab_imp)
avg_imp <- apply(split_imp,1,mean)
num_keep <- ceiling(length(avg_imp)*prop_top)
top_rats <- names(avg_imp)[pick_top_n(avg_imp,num_keep)]
writeLines(top_rats,paste(push_iter,"_ratios.txt",sep=""))
for(i in mod_list){file.remove(paste(path1,i,path2,sep=""))}
