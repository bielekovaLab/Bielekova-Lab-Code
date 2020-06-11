my_filename <- function(string){
  the_date <- str_replace_all(as.character(Sys.Date()),"-","")
  split_string <- unlist(strsplit(string,".",fixed=TRUE))
  n <- length(split_string)
  split_string[n-1] <- paste(split_string[n-1],the_date,sep="_")
  out <- paste(split_string,collapse = ".")
  return(out)
}
clean <- function(dat,...){
  # dat <- clinical
  r_miss <- apply(dat,1,function(x){!all(is.na(x))})
  c_miss <- apply(dat,2,function(x){!all(is.na(x))})
  dat <- dat[r_miss,c_miss]
  dat <- clean_names(dat)
  if(any(names(dat) %in% top_order)){
    ind <-  names(dat) %in% top_order
    new_order <- top_order[sort(match(names(dat)[ind],top_order))]
    new_order <- c(new_order,names(dat)[!ind])
    dat <- dat[,new_order]
  }
  if(all(c("diagnosis_main","diagnosis_category") %in% names(dat))){
    dat <- mutate(dat,diagnosis = clean_diagnosis(diagnosis_main,diagnosis_category)) %>% 
      select(-diagnosis_main,-diagnosis_category)
  }
  # if(all(c("patientcode","date","protocol1","neurex_total","neurex_total_thadj","neurex_rating_completed") %in% names(dat))){
  #   dat <- clean_neurex(dat)
  # }
  # if(all(c("patientcode","date","protocol1","combiwise_score","combiwise_score_thadj","combiwise_validation") %in% names(dat))){
  #   dat <- clean_combiwise(dat)
  # }
  return(dat)
}
cohort_order <- c("discovery","validation","DAC-HYP","denver","Rivitalise MedImmune",
                  "June 2015","Nov 2015","May 2016","CHI 1.3k MAY '17","CHI October 2018")
clean_names <- function(dat,...){
  # dat <- trials
  d_names <- names(dat)
  d_names <- sapply(strsplit(d_names,"::"),last)
  d_names <- tolower(d_names)
  d_names <- str_replace(d_names," ","_")
  d_names <- str_replace(d_names,"-","_")
  names(dat) <- d_names
  return(dat)
}
top_order <- c("sampleid","patientcode","patientid","age","date","lpdate","visit","periodnumber",
               "protocol1","therapy","other_dmts")
clean_diagnosis <- function(diagnosis_main, diagnosis_category){
  # diagnosis_main <- diagnosis$diagnosis_main
  # diagnosis_category <- diagnosis$diagnosis_category
  diagnosis_main <- sapply(strsplit(diagnosis_main," "),first)
  diagnosis_main <- ifelse(diagnosis_main=="Healthy","HD",diagnosis_main)
  diagnosis <- ifelse(is.na(diagnosis_main),NA,ifelse(diagnosis_main=="MS",diagnosis_category,diagnosis_main))
  return(diagnosis)
}

# clean_combiwise, only works for raw combiwise and therapy adjusted combiwise, if both are available
# and patientcode, date, and protocol1 are available
clean_combiwise <- function(dat){
  # dat <- clinical
  combi_dat <- dat %>% 
    filter(is.na(combiwise_validation)) %>% 
    select(patientcode,date,protocol1,combiwise_score,combiwise_score_thadj)
  nocombi_dat <- dat %>% 
    select(-combiwise_score,-combiwise_score_thadj,-combiwise_validation)
  out_dat <- merge(nocombi_dat,combi_dat,by=c("patientcode","date","protocol1"),all.x=TRUE)
  return(out_dat)
}

# clean neurex, only works for neurex_total and neurex_total_thadj, if both are available
# and patientcode, date, and protocol1 are available
clean_neurex <- function(dat){
  # dat <- clinical
  neurex_dat <- dat %>% 
    filter(!is.na(neurex_rating_completed)) %>% 
    filter(neurex_rating_completed==1) %>% 
    select(patientcode,date,protocol1,neurex_total,neurex_total_thadj)
  noneurex_dat <- dat %>% 
    select(-neurex_total,-neurex_total_thadj,-neurex_rating_completed)
  out_dat <- merge(noneurex_dat,neurex_dat,by=c("patientcode","date","protocol1"),all.x=TRUE)
  return(out_dat)
}


match_dates <- function(dat1,dat2,date_name,g_ind,t_ind){
  # dat1 <- old_dat %>%
  #   select(patient,lpdate,sampleid)
  # dat2 <- clinical
  # date_name <- "clinical_date"
  # g_ind <- 1
  # t_ind <- 2
  p_ind <- which(dat1[[g_ind]] %in% dat2[[g_ind]])
  d1 <- split(dat1,dat1[[g_ind]])
  d2 <- split(dat2,dat2[[g_ind]])
  g1 <- names(d1)
  g2 <- names(d2)
  for(i in 1:length(g1)){
    if(length(which(g2==g1[i]))==0){
      d1[[i]][date_name] <- rep(NA,dim(d1[[i]])[1])
    }
    if(length(which(g2==g1[i]))!=0){
      sd1 <- d1[[g1[i]]]
      sd2 <- d2[[which(g2==g1[i])]]
      date1 <- sd1[[t_ind]]
      date2 <- sd2[[t_ind]]
      to_merge <- rep(NA,length=length(date1))
      for(j in 1:length(date1)){
        diffs <- as.numeric(abs(date1[j] - date2))
        to_merge[j] <- which.min(diffs)
      }
      d1[[i]][date_name] <- as.character(date2[to_merge])
    }
  }
  out <- foreach(i=1:length(d1),.combine='rbind') %do% tbl_df(d1[[i]])
  return(out)
}
ratio_append <- function(data,ratio_list,expo=FALSE,adj=FALSE){
  # data <- old_untreated
  # ratio_list <- var_lists[[4]]
  splitlist <- strsplit(ratio_list,"[/]")
  v1 <- sapply(splitlist,first)
  v2 <- sapply(splitlist,last)
  p <- length(ratio_list)
  
  if(adj==TRUE){for(i in 1:p){
    data[,ratio_list[i]] <- data[[v1[i]]] - data[[v2[i]]]
  }}
  if(expo==FALSE){for(i in 1:p){
    data[,ratio_list[i]] <- data[[v1[i]]]/data[[v2[i]]]
  }}
  if(expo==TRUE){for(i in 1:p){
    data[,ratio_list[i]] <- exp(data[[v1[i]]])/exp(data[[v2[i]]])
  }}
  return(data)
}
