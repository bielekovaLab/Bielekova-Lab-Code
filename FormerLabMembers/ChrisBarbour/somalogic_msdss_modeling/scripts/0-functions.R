library(readr, quietly = TRUE)
library(readxl, quietly = TRUE)
# library(plyr)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(purrr, quietly = TRUE)
library(bindrcpp, quietly = TRUE)
library(data.table, quietly = TRUE)
library(foreach, quietly = TRUE)
library(stringr, quietly = TRUE)
library(ranger)
library(parallel)
library(doParallel)
library(batch)
library(car)
library(gbm)


detectBatchCPUs <- function() { 
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
  if (is.na(ncores)) { 
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
  } 
  if (is.na(ncores)) { 
    return(4) # for helix
  } 
  return(ncores) 
}
adjust_function <- function(y_char,data=new,alpha_val=0.05,...){
  # y_char <- "SL000001/SL000002"
  # y_char <- "SL000001"
  # data <- new
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

# adjust_function <- function(y_char,data=new,alpha_val=0.05,...){
#   # y_char <- "SL000001/SL000002"
#   # y_char <- "SL000001"
#   # data <- new
#   if(nchar(y_char)==17){data <- ratio_append(data,y_char,expo = FALSE)}
#   y<-log(data[[y_char]])
#   # mod <- lm(y~age*gender,data)
#   if(TRUE){
#     mod <- lm(y~age+gender,data)
#     pvals <- Anova(mod,type=2)$`Pr(>F)`[1:2]
#     if(any(pvals>alpha_val)){
#       ind <- which(pvals>alpha_val)
#       if(length(ind)==2){
#         mod <- lm(y~1)
#       }
#       if(length(ind)==1){
#         if(ind==1){
#           mod <- lm(y~gender,data)
#           if(anova(mod)$`Pr(>F)`[1]>alpha_val){
#             mod <- lm(y~1)
#           }
#         }
#         if(ind==2){
#           mod <- lm(y~age,data)
#           if(anova(mod)$`Pr(>F)`[1]>alpha_val){
#             mod <- lm(y~1)
#           }
#         }
#       }
#     }
#   }
#   # y_adjust <- residuals(mod)
#   # return(coef(mod))
#   return(mod)
# }
char_2_fac <- function(x){
  if(x != "character"){return(x)}
  if(x == "character"){return(as.factor(x))}
}
create_rat <- function(rat,dat=soma,expo=FALSE,adj=FALSE){
  if(nchar(rat)==17){
    ind_mark <- split_rat(rat)
    if(adj==TRUE){return(dat[[ind_mark[1]]] - dat[[ind_mark[2]]])}
    if(expo==FALSE & adj==FALSE){return(dat[[ind_mark[1]]]/dat[[ind_mark[2]]])}
    if(expo==TRUE){return(exp(dat[[ind_mark[1]]])/exp(dat[[ind_mark[2]]]))}
  }
  if(nchar(rat) != 17){
    return(dat[[rat]])
  }
}
split_rat <- function(rat){
  split_list <- unlist(strsplit(rat,"[/]"))
  v1 <- first(split_list)
  v2 <- last(split_list)
  return(c(v1,v2))
}
quant_bin <- function(x,quants=c(0,0.5,1),...){
  quantiles <- quantile(x,quants)
  quantiles[1] <- quantiles[1]-0.001
  out <- cut(x,breaks=quantiles,...)
  return(out)
}
my_cor <- function(x,y,...){
  obj <- cor.test(x,y,...)
  val <- obj$estimate
  pval <- obj$p.value
  out <- paste(round(val,3)," (",format.pval(pval,eps=0.001,digits=4),")",sep="")
}

remove_dups <- function(dat){
  dat <- dat[!duplicated(dat),]
  return(dat)
}
grab_soma <- function(mark_frame,path="./data/translation file.xlsx",platform="1.1K"){
  # path <- "./data/translation file.xlsx"
  # mark_vec <- marker_overlap[,1:2]
  # platform <- "1.1K"
  # mark_frame<- df_frame1
  tomatch <- read_excel(path,col_names = FALSE)
  if(platform == "1.3K"){
    tomatch <- tomatch[,2:4]
  }
  if(platform == "1.1K"){
    tomatch <- tomatch[,11:13]
  }
  names(tomatch) <- c("soma","gene","protein")
  gene_mean <- mean(mark_frame[[1]] %in% tomatch[[2]])
  protein_mean <- mean(mark_frame[[2]] %in% tomatch[[3]])
  if(gene_mean != 1|protein_mean!=1){stop("Something not quite right.... check other platform!")}
  d_frame <- data.frame(gene=mark_frame[[1]],protein=mark_frame[[2]])
  matched <- tbl_df(merge(d_frame,tomatch,by=c("gene","protein"),sort = FALSE))
  return(matched[["soma"]])
}
grab_soma_rat <- function(ratio_frame,path="./data/translation file.xlsx",platform="1.3K"){
  # path <- "./data/translation file.xlsx"
  # ratio_frame <- ratio_overlap[,1:2]
  # platform <- "1.3K"
  tomatch <- read_excel(path,col_names = FALSE)
  if(platform == "1.3K"){
    tomatch <- tomatch[,2:4]
  }
  if(platform == "1.1K"){
    tomatch <- tomatch[,11:13]
  } 
  names(tomatch) <- c("soma","gene","protein")
  tomatch <- arrange(tomatch,soma)
  genes <- tomatch[["gene"]]
  proteins <- tomatch[["protein"]]
  soma <- tomatch[["soma"]]
  n <- length(genes)
  gene_ratios <- soma_ratios <- protein_ratios <- choose(n,2)
  iter <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      soma_ratios[iter] <- paste(soma[i],soma[j],sep="/")
      gene_ratios[iter] <- paste(genes[i],genes[j],sep="/")
      protein_ratios[iter] <- paste(proteins[i],proteins[j],sep="/")
      iter <- iter + 1
    }
  }
  tomatch1 <- tibble(soma=soma_ratios,gene=gene_ratios,protein=protein_ratios)
  gene_ratios <- soma_ratios <- protein_ratios <- choose(n,2)
  iter <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      soma_ratios[iter] <- paste(soma[j],soma[i],sep="/")
      gene_ratios[iter] <- paste(genes[j],genes[i],sep="/")
      protein_ratios[iter] <- paste(proteins[j],proteins[i],sep="/")
      iter <- iter + 1
    }
  }
  tomatch2 <- tibble(soma=soma_ratios,gene=gene_ratios,protein=protein_ratios)
  tomatch <- rbind(tomatch1,tomatch2)
  gene_mean <- mean(ratio_frame[[1]] %in% tomatch[[2]])
  protein_mean <- mean(ratio_frame[[2]] %in% tomatch[[3]])
  if(gene_mean != 1|protein_mean!=1){stop("Something not quite right.... check other platform!")}
  d_frame <- data.frame(gene=ratio_frame[[1]],protein=ratio_frame[[2]])
  matched <- tbl_df(merge(d_frame,tomatch,by=c("gene","protein"),sort = FALSE))
  return(matched[["soma"]])
}
clean_rat <- function(rat_vec,path="./data/translation file.xlsx",platform="1.3K"){
  # path <- "./data/translation file.xlsx"
  # rat_vec <- var_lists[[1]]
  # platform <- "1.1K"
  tomatch <- read_excel(path,col_names = FALSE)
  if(platform == "1.3K"){
    tomatch <- tomatch[,2:4]
  }
  if(platform == "1.1K"){
    tomatch <- tomatch[,11:13]
  }
  names(tomatch) <- c("soma","gene","protein")
  rat_dat <- data.frame(ratio=as.character(rat_vec)) %>% 
    tbl_df() %>% 
    mutate(ratio = toupper(as.character(ratio))) %>% 
    mutate(v1 = sapply(strsplit(ratio,"/"),dplyr::first)) %>% 
    merge(tomatch,by.x="v1",by.y="soma",sort=FALSE) %>% 
    rename(gene1 = gene,protein1=protein) %>% 
    mutate(v2 = sapply(strsplit(ratio,"/"),dplyr::last)) %>% 
    merge(tomatch,by.x="v2",by.y="soma",sort=FALSE) %>% 
    rename(gene2 = gene,protein2=protein) %>% 
    dplyr::select(ratio,v1,gene1,protein1,v2,gene2,protein2) %>% 
    arrange(ratio)
  return(rat_dat)
}
clean_marker <- function(marker_vec,path="./data/translation file.xlsx",platform="1.3K"){
  # path <- "./data/translation file.xlsx"
  # marker_vec <- names(which(apply(model_dat[,msdss_last_rats],2,function(x){sum(is.na(x))}) !=0))[-1]
  # platform <- "1.1K"
  tomatch <- read_excel(path,col_names = FALSE)
  if(platform == "1.3K"){
    tomatch <- tomatch[,2:4]
  }
  if(platform == "1.1K"){
    tomatch <- tomatch[,11:13]
  }
  names(tomatch) <- c("soma","gene","protein")
  marker_dat <- data.frame(marker=as.character(marker_vec)) %>% 
    tbl_df() %>% 
    mutate(marker = as.character(marker)) %>% 
    # mutate(v1 = sapply(strsplit(ratio,"/"),dplyr::first)) %>% 
    merge(tomatch,by.x="marker",by.y="soma",sort=FALSE) %>% 
    rename(gene = gene,protein=protein) %>% 
    # mutate(v2 = sapply(strsplit(ratio,"/"),dplyr::last)) %>% 
    # merge(tomatch,by.x="v2",by.y="soma") %>% 
    # rename(gene2 = gene,protein2=protein) %>% 
    dplyr::select(marker,gene,protein)
  return(marker_dat)
}


rat_to_gene <- function(rat_vec,path="./data/translation file.xlsx",platform="1.3K"){
  tomatch <- read_excel(path,col_names = FALSE)
  if(platform == "1.3K"){
    tomatch <- tomatch[,2:4]
  }
  if(platform == "1.1K"){
    tomatch <- tomatch[,11:13]
  }
  names(tomatch) <- c("soma","gene","protein")
  rat_dat <- data.frame(ratio=as.character(rat_vec)) %>% 
    tbl_df() %>% 
    mutate(ratio = toupper(as.character(ratio))) %>% 
    mutate(v1 = sapply(strsplit(ratio,"/"),dplyr::first)) %>% 
    merge(tomatch,by.x="v1",by.y="soma",sort=FALSE) %>% 
    rename(gene1 = gene,protein1=protein) %>% 
    mutate(v2 = sapply(strsplit(ratio,"/"),dplyr::last)) %>% 
    merge(tomatch,by.x="v2",by.y="soma",sort=FALSE) %>% 
    rename(gene2 = gene,protein2=protein) %>% 
    dplyr::select(ratio,v1,gene1,protein1,v2,gene2,protein2)
  out <- with(rat_dat,paste(gene1,"/",gene2,sep=""))
  return(out)
}
marker_to_gene <- function(marker_vec,path="./data/translation file.xlsx",platform="1.3K"){
  # path <- "./data/translation file.xlsx"
  # marker_vec <- names(which(apply(model_dat[,msdss_last_rats],2,function(x){sum(is.na(x))}) !=0))[-1]
  # platform <- "1.1K"
  tomatch <- read_excel(path,col_names = FALSE)
  if(platform == "1.3K"){
    tomatch <- tomatch[,2:4]
  }
  if(platform == "1.1K"){
    tomatch <- tomatch[,11:13]
  }
  names(tomatch) <- c("soma","gene","protein")
  marker_dat <- data.frame(marker=as.character(marker_vec)) %>% 
    tbl_df() %>% 
    mutate(marker = as.character(marker)) %>% 
    merge(tomatch,by.x="marker",by.y="soma",sort=FALSE) %>% 
    rename(gene = gene,protein=protein) %>% 
    dplyr::select(marker,gene,protein)
  return(marker_dat$gene)
}

'%notin%' <- function(x,y)!('%in%'(x,y))
